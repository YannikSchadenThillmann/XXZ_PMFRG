
function fng(Par) #filenamegenerator
    Name = Par.System.Name
    (;N,theta) = Par.NumericalParams
    (;NLen,couplings) = Par.System
    
    coupstring = ""
    for i in eachindex(couplings)
        if couplings[i] != 0.0
            coupstring = coupstring * "J$(i-1)=$(couplings[i])_"
        end
    end

    return "PMFRG_XXZ_Tflow_$(Name)_$(coupstring)N=$(N)_NLen=$(NLen)_theta=$(string(round(theta/pi,digits = 2)))pi.h5"
end

function stround(number::Float64, digits::Int64)
    res = string(round(number, digits = digits))
    digits = digits + length(split(res,".")[1]) + 1
    reps = digits - length(res)
    if reps > 0
        return res*repeat("0",reps)
    else
        return res
    end
end

function stround(number::String, digits::Int64)
    res = string(round(parse(Float64,number), digits = digits))
    digits = digits + length(split(res,".")[1]) + 1
    reps = digits - length(res)
    if reps > 0
        return res*repeat("0",reps)
    else
        return res
    end
end

function countdigits(n::Float64) # counting significant digits after floating point
    return Int64(ceil(abs(log10(n))))
end

function h5write(filename, name::AbstractString, data; pv...)
    file = h5open(filename, "cw")
    try
        write(file, name, data)
    finally
        close(file)
    end
end

function uniquename(path::String, print::Bool)
    posfolder = findlast('/',path)
    posext = findlast('.',path)
    folder = path[1:posfolder]
    filename = path[posfolder+1:posext-1]
    ext = path[posext:end]

    dir = readdir(folder)
    n = 1
    while filename*ext in dir || "Checkpoints_"*filename in dir
        if print == true
            if filename*ext in dir
                println("File ", filename * ext," already exists in ", folder)
            elseif "Checkpoints_"*filename in dir
                println("Folder ", "Checkpoints_"*filename," already exists in ", folder)
            end
        end

        if findlast("_copy",filename) == nothing
            filename = filename*"_copy$(n)"
        else
            n += 1
            filename = filename[1:end-1]*string(n)
        end
    end

    return (folder, filename, ext)
end
uniquename(::Nothing, args...) = nothing

function saveMainOutput(filename::String, saved_values, Par, Katanin)
    println("Saving Main output to ", filename)

    file = h5open(filename, "cw")

    ObsArr = StructArray(saved_values.saveval)
    fields = fieldnames(Observables)
    for f in fields
        Arr = convertToArray(getproperty(ObsArr,f))
        write(file,"State/"*string(f),Arr)
    end

    fields = fieldnames(NumericalParamsXXZ)
    for f in fields
        write(file,"Par/NumericalParams/"*string(f),getproperty(Par.NumericalParams,f))
    end

    fields = (:Name, :NLen, :couplings, :Npairs)
    for f in fields
        write(file,"Par/System/"*string(f),getproperty(Par.System,f))
    end
    
    write(file,"State/Tvals",saved_values.t)
    write(file,"Par/Katanin",Katanin)
    close(file)
end
saveMainOutput(::Nothing,args...) = nothing

function convertToArray(VecOfArray::AbstractVector{VT}) where {N,VT <: AbstractArray{T,N} where T}
    cat(VecOfArray...,dims = N+1)
end

function getObservables(::Type{Observables},State::ArrayPartition,T,Par)
    (;N) = Par.NumericalParams
    f_int,gammax,gammaz,Vxxxx,Vzzzz,Vxxyy,Vxxzz,Vzzxx,Vxyxy,Vxzxz,Vzxzx = State.x
    chixx = getChixx(State,T,Par,N)
    chizz = getChizz(State,T,Par,N)
    MaxVxxxx = maximum(abs,Vxxxx,dims = (2,3,4,5))[:,1,1,1]
    MaxVzzzz = maximum(abs,Vzzzz,dims = (2,3,4,5))[:,1,1,1]
    MaxVxxyy = maximum(abs,Vxxyy,dims = (2,3,4,5))[:,1,1,1]
    MaxVxxzz = maximum(abs,Vxxzz,dims = (2,3,4,5))[:,1,1,1]
    MaxVzzxx = maximum(abs,Vzzxx,dims = (2,3,4,5))[:,1,1,1]
    MaxVxyxy = maximum(abs,Vxyxy,dims = (2,3,4,5))[:,1,1,1]
    MaxVxzxz = maximum(abs,Vxzxz,dims = (2,3,4,5))[:,1,1,1]
    MaxVzxzx = maximum(abs,Vzxzx,dims = (2,3,4,5))[:,1,1,1]
    return Observables(chixx,chizz,copy(gammax),copy(gammaz),copy(f_int),MaxVxxxx,MaxVzzzz,MaxVxxyy,MaxVxxzz,MaxVzzxx,MaxVxyxy,MaxVxzxz,MaxVzxzx) # make sure to allocate new memory each time this function is called
end

function saveCheckpoint(mainfile::Tuple{String, String, String}, saved_values, Par, i::Int64, dtmin)
    (folder, file, ext) = mainfile

    CPfolder = folder*"Checkpoints_$(file)/"
    !("Checkpoints_$(file)" in readdir(folder)) && (mkpath(CPfolder); println("Creating path: ", CPfolder))

    file = CPfolder*"Checkpoint_$(i)_$(file)_T=$(stround(t_to_T(saved_values.t[end]), countdigits(dtmin)+1))"*ext

    println("Saving Checkpoint to ", file)

    file = h5open(file, "cw")
    ObsArr = StructArray(saved_values.saveval)
    fields = fieldnames(Observables)
    for f in fields
        Arr = convertToArray(getproperty(ObsArr,f))
        write(file,"State/"*string(f),Arr)
    end

    fields = fieldnames(NumericalParamsXXZ)
    for f in fields
        write(file,"Par/NumericalParams/"*string(f),getproperty(Par.NumericalParams,f))
    end

    fields = (:Name, :NLen, :couplings, :Npairs)
    for f in fields
        write(file,"Par/System/"*string(f),getproperty(Par.System,f))
    end
    
    write(file,"State/Tvals",t_to_T.(saved_values.t))
    close(file)

    if i > 3
        oldfile = filter(x->occursin("Checkpoint_$(i-3)",x), readdir(CPfolder))[1]
        rm(CPfolder*oldfile)
    end
end
saveCheckpoint(::Nothing,args...) = nothing
