
function getChannel(Buffs::AbstractVector{<:T}) where {T}
    BufferChannel = Channel{T}(length(Buffs))
    for buff in Buffs
        put!(BufferChannel, buff)
    end
    return BufferChannel
end

function AllocateSetup(Par::OneLoopParamsXXZ)
    println("XXZ PMFRG, One Loop, Tflow, Jx/Jz = tan(θ="*string(round(Par.NumericalParams.theta/pi, digits = 4))*"π) ≈ "*string(round(tan(Par.NumericalParams.theta), digits = 2)))
    ##Allocate Memory:
    X = BubbleTypeXXZ(Par)

    VertexBuffers =
    getChannel([VertexBufferType(typeof(Par.NumericalParams.T_min), Par.System.Npairs) for _ = 1:Threads.nthreads()])

    return (;X,Par,VertexBuffers)
end

function InitializeState(Par)
    (;N,Ngamma,theta) = Par.NumericalParams
    VDims = getVDimsXXZ(Par)
    (;couplings,NUnique) = Par.System

    floattype = _getFloatType(Par)
    
    State = ArrayPartition( #Allocate Memory:
        zeros(floattype,NUnique), # f_int 
        zeros(floattype,NUnique,Ngamma), # gamma.x
        zeros(floattype,NUnique,Ngamma), # gamma.z
        zeros(floattype,VDims), # V_xxxx
        zeros(floattype,VDims), # V_zzzz
        zeros(floattype,VDims), # V_xxyy
        zeros(floattype,VDims), # V_xxzz
        zeros(floattype,VDims), # V_zzxx
        zeros(floattype,VDims), # V_xyxy
        zeros(floattype,VDims), # V_xzxz
        zeros(floattype,VDims)  # V_zxzx
    )

    Γxyxy = State.x[9]
    Γxzxz = State.x[10]
    Γzxzx = State.x[11]
    setToBareVertex!(Γxyxy,Γxzxz,Γzxzx,couplings,theta)
    return State
end

function setToBareVertex!(Γxyxy::AbstractArray{T,4},Γxzxz::AbstractArray{T,4},Γzxzx::AbstractArray{T,4},couplings::AbstractVector, theta::Float64) where T
    Jx = sin(theta)
    Jz = cos(theta)
    for Rj in axes(Γxyxy,1)
        Γxyxy[Rj,:,:,:] .= -Jz.*couplings[Rj]
        Γxzxz[Rj,:,:,:] .= -Jx.*couplings[Rj]
        Γzxzx[Rj,:,:,:] .= -Jx.*couplings[Rj]
    end
    return Γxyxy, Γxzxz, Γzxzx
end

_getFloatType(Par) = typeof(Par.NumericalParams.T_min)

SolveFRGXXZ(Par;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)

function launchPMFRG!(State, setup, Deriv!::Function;
    method = DP5(),
    ObsSaveat = nothing,
    ObservableType = Observables,
    mainfile = nothing,
    checkpointinterval = 1,
    dtmin = 1e-4,
    Katanin = true,
    kwargs...
    )

    Par = setup[2]
    (;T_max,T_min,accuracy) = Par.NumericalParams

    if Katanin
        println("Katanin is on")
    else
        println("Katanin is off")
    end

    # unique name for file and checkpoint folder
    mainfile = uniquename(mainfile, false)

    # checkpoints
    i = 0
    function Checkpoint(state,t,integrator)
        i+=1
        if i%checkpointinterval == 0
            saveCheckpoint(mainfile, saved_values, Par, i, dtmin)
        end
        println("t = $t")
    end

    # define Callbacks (points in flow, where quantities are computed)
    ObsSaveat = gettMesh(ObsSaveat,T_min,T_max)
    save_func(State,t,integrator) = getObservables(ObservableType,State,t_to_T(t),Par)
    saved_values = SavedValues(eltype(State),ObservableType)
    saveCB = SavingCallback(save_func, saved_values, save_everystep =false, saveat = ObsSaveat, tdir=-1)
    checkpointCB = FunctionCallingCallback(Checkpoint,tdir=-1,func_start = false)

    # ODE solver
    t0 = T_to_t(T_max)
    tend = get_t_min(T_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!, Katanin)
    problem = ODEProblem(Deriv_subst!,State,(t0,tend),setup)
    sol = solve(problem, method, reltol = accuracy, abstol = accuracy, save_everystep = false, callback=CallbackSet(saveCB,checkpointCB), dt=T_to_t(0.2*T_max), dtmin=dtmin; kwargs...)

    # save quantities
    saved_values.t .= t_to_T.(saved_values.t)
    if mainfile !== nothing
        saveMainOutput(join(mainfile), saved_values, Par, Katanin)
    end
    return sol
end

t_to_T(t) = exp(t)
T_to_t(T) = log(T)

function get_t_min(T)
    T < exp(-30) && @warn "T_min too small! Set to exp(-30) instead."
    max(T_to_t(T),-30.)
end

# instead of T, formulate differential equation in log(T)
function generateSubstituteDeriv(getDeriv!::Function, Katanin)
    
    function DerivSubs!(Deriv,State,par,t)
        T = t_to_T(t)
        a = getDeriv!(Deriv,State,par,T, Katanin)
        Deriv .*= T
        a
    end
end

# determine points in flow, at which Callbacks are performed (where quantities are saved)
function getTempMesh(Saveat::Nothing,T_min,T_max)
    dense_range = collect(LinRange(T_min,2.,400))
    medium_range = collect(LinRange(2.,10.,100))
    sparse_range = collect(LinRange(10.,T_max,30))
    ObsSaveat = unique!(append!(dense_range,medium_range,sparse_range))
    return ObsSaveat
end

function getTempMesh(Saveat::Vector{Float64},T_min,T_max)
    return unique(push!(Saveat,T_max)) # make sure that there is at least one element at beginning of code
end

gettMesh(Saveat,T_min,T_max) = T_to_t.(getTempMesh(Saveat,T_min,T_max))