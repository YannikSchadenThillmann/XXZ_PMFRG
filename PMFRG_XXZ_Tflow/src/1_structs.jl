
# auxiliary functions
function test(a::String)
    println(a)
end

strd(x,n) = string(round(x,digits=n))

setZero!(a::AbstractArray{T,N}) where {T,N} = fill!(a,zero(T))

function setZero!(PartArr::ArrayPartition)
    for arr in PartArr.x
        fill!(arr,0.)
    end
end

"""Recursively sets structure to zero"""
function setZero!(a::T) where T
    for f in fieldnames(T)
        setZero!(getfield(a,f))
    end
    return a
end


# structs for parameters, vertices and states
struct NumericalParamsXXZ
    theta::Float64 # Jx = sin(theta) and Jz = cos(theta)
    N::Int
    Ngamma::Int

    accuracy::Float64
    T_min::Float64
    T_max::Float64

    lenIntw::Int
    lenIntw_acc::Int
    np_vec::Array{Int,1}
    np_vec_gamma::Array{Int,1}
end

function NumericalParamsXXZ(;
    theta::Float64 = 0.25*pi, # isotropic case, if theta is not specified
    N::Integer = 24,
    Ngamma::Integer = N, #Number of gamma frequencies

    accuracy = 1e-6,
    T_min = exp(-10.),
    T_max = exp(10.),

    lenIntw::Int = N,
    lenIntw_acc::Int = 2*maximum((N,Ngamma,lenIntw)),
    np_vec::Array{Int,1} = collect(0:N-1),
    np_vec_gamma::Array{Int,1} = collect(0:Ngamma-1)
    )

    return NumericalParamsXXZ(
        theta,
        N,
        Ngamma,

        accuracy,
        T_min,
        T_max,

        lenIntw,
        lenIntw_acc,
        np_vec,
        np_vec_gamma,
    )
end

struct SelfEnergyType{T}
    x::Array{T,2}
    z::Array{T,2}
end

function SelfEnergyType(SEDims::Tuple)
    return SelfEnergyType(
        zeros(Float64,SEDims), # gamma_x
        zeros(Float64,SEDims) # gamma_z
    )
end
getSEDims(Par) = (Par.System.NUnique,Par.NumericalParams.Ngamma)
SelfEnergyType(Par) = SelfEnergyType(getSEDims(Par))

struct VertexTypeXXZ{T}
    xxxx::Array{T,4}
    zzzz::Array{T,4}

    xxyy::Array{T,4}
    xxzz::Array{T,4}
    zzxx::Array{T,4}

    xyxy::Array{T,4}
    xzxz::Array{T,4}
    zxzx::Array{T,4}
end

function VertexTypeXXZ(VDims::Tuple)
    return VertexTypeXXZ(
        zeros(VDims), # Gamma_xxxx
        zeros(VDims), # ...

        zeros(VDims),
        zeros(VDims),
        zeros(VDims),

        zeros(VDims),
        zeros(VDims),
        zeros(VDims)
    )
end
getVDimsXXZ(Par) = (Par.System.Npairs,Par.NumericalParams.N,Par.NumericalParams.N,Par.NumericalParams.N)
VertexTypeXXZ(Par) = VertexTypeXXZ(getVDimsXXZ(Par))

struct StateTypeXXZ{T}
    f_int::Array{T,1}
    γ::SelfEnergyType{T}
    Γ::VertexTypeXXZ{T}
end

function StateTypeXXZ(NUnique::Int,VDims::Tuple,SEDims::Int64,type = Float64::Type)
    return StateTypeXXZ(
        zeros(type,NUnique), # f int
        SelfEnergyType(SEDims), # Self Energy Dimensions for gamma.x and gamma.z
        VertexTypeXXZ(VDims)
    )
end

StateTypeXXZ(Par) = StateTypeXXZ(Par.System.NUnique,Par.NumericalParams.Ngamma,getVDimsXXZ(Par),_getFloatType(Par))
StateTypeXXZ(f_int,γx,γz,Γxxxx,Γzzzz,Γxxyy,Γxxzz,Γzzxx,Γxyxy,Γxzxz,Γzxzx) = StateTypeXXZ(f_int,SelfEnergyType(γx,γz),VertexTypeXXZ(Γxxxx,Γzzzz,Γxxyy,Γxxzz,Γzzxx,Γxyxy,Γxzxz,Γzxzx)) 

RecursiveArrayTools.ArrayPartition(x) = ArrayPartition(x.f_int,x.γ.x,x.γ.z,x.Γ.xxxx,x.Γ.zzzz,x.Γ.xxyy,x.Γ.xxzz,x.Γ.zzxx,x.Γ.xyxy,x.Γ.xzxz,x.Γ.zxzx)
StateTypeXXZ(Arr::ArrayPartition) = StateTypeXXZ(Arr.x...)

struct BubbleTypeXXZ{T}
    # X bubbles
    xxxx::Array{T,4}
    zzzz::Array{T,4}

    xxyy::Array{T,4}
    xxzz::Array{T,4}
    zzxx::Array{T,4}

    xyxy::Array{T,4}
    xzxz::Array{T,4}
    zxzx::Array{T,4}

    # XTilde bubbles
    Txxxx::Array{T,4}
    Tzzzz::Array{T,4}

    Txxyy::Array{T,4}  
    Txxzz::Array{T,4}
    Tzzxx::Array{T,4}

    Txyxy::Array{T,4}
    Txzxz::Array{T,4}
    Tzxzx::Array{T,4}

    Txyyx::Array{T,4} # additional 3 vertices for XTilde
    Txzzx::Array{T,4}
    Tzxxz::Array{T,4}
end

function BubbleTypeXXZ(VDims::Tuple,type = Float64)
    return BubbleTypeXXZ(
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),

        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims),
        zeros(type,VDims)
    )
end
BubbleTypeXXZ(Par) = BubbleTypeXXZ(getVDimsXXZ(Par))

struct OptionParamsXXZ
    usesymmetry::Bool
    MinimalOutput::Bool
end
OptionParamsXXZ(;usesymmetry::Bool = false, MinimalOutput::Bool = false,kwargs...) = OptionParamsXXZ(usesymmetry,MinimalOutput)

struct OneLoopParamsXXZ{S,N,O}
    System::S
    NumericalParams::N
    Options::O
end

ParamsXXZ(System;kwargs...) = OneLoopParamsXXZ(System,NumericalParamsXXZ(;kwargs...),OptionParamsXXZ(;kwargs...))

struct OneLoopWorkspaceXXZ{T,ParType}
    State::StateTypeXXZ{T}  #Stores the current state
    Deriv::StateTypeXXZ{T}  #Stores the derivative
    X::BubbleTypeXXZ{T} #Stores the bubble function X and XTilde
    Par::ParType # Params
end

function OneLoopWorkspaceXXZ(Deriv,State,X,Par)
    setZero!(Deriv)
    setZero!(X)
    return OneLoopWorkspaceXXZ(
        StateTypeXXZ(State.x...),
        StateTypeXXZ(Deriv.x...),
        X,
        Par
    )
end

struct Observables{T}
    chixx::Matrix{T}
    chizz::Matrix{T}
    gammax::Matrix{T}
    gammaz::Matrix{T}
    f_int::Vector{T}
    maxVxxxx::Vector{T}
    maxVzzzz::Vector{T}
    maxVxxyy::Vector{T}
    maxVxxzz::Vector{T}
    maxVzzxx::Vector{T}
    maxVxyxy::Vector{T}
    maxVxzxz::Vector{T}
    maxVzxzx::Vector{T}
end

struct VertexBufferType{T}
    V_xxxx_21::Vector{T}
    V_zzzz_21::Vector{T}
    V_xxyy_21::Vector{T}
    V_xxzz_21::Vector{T}
    V_zzxx_21::Vector{T}
    V_xyxy_21::Vector{T}
    V_xzxz_21::Vector{T}
    V_zxzx_21::Vector{T}
    V_xyxy_12::Vector{T}
    V_xzxz_12::Vector{T}
    V_zxzx_12::Vector{T}
    V_xxxx_34::Vector{T}
    V_zzzz_34::Vector{T}
    V_xxyy_34::Vector{T}
    V_xxzz_34::Vector{T}
    V_zzxx_34::Vector{T}
    V_xyxy_34::Vector{T}
    V_xzxz_34::Vector{T}
    V_zxzx_34::Vector{T}
    V_xyxy_43::Vector{T}
    V_xzxz_43::Vector{T}
    V_zxzx_43::Vector{T}
end
VertexBufferType(type, Npairs) = VertexBufferType((zeros(type, Npairs) for _ = 1:22)...)