
# cutoff function
θ(T) = 1/sqrt(T)

# 1 particle Propagators
function get_w(nw)
    return (2*nw+1)*pi
end

function iG0(nw,T)
    w = get_w(nw)
    θ(T)/(w)
end

function iG_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer)
    w = get_w(nw)
    return 1/(w/θ(T) + gamma_(gamma,x,nw))
end

function iS_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer)
    w = get_w(nw)
    return -θ(T) *w/2 *iG_(gamma,x,T,nw)^2
end

function iSKat_(gamma::AbstractArray, Dgamma::AbstractArray, x::Integer, T::Real, nw::Integer, Katanin::Bool)
    w = get_w(nw)
    return -(θ(T) *w/2 + Float64(Katanin)* gamma_(Dgamma,x,nw)) *iG_(gamma,x,T,nw)^2
end

function gamma_(gamma::AbstractArray, x::Integer, nw::Integer)
    Ngamma = size(gamma,2)
    s = 1
    if nw<0
        nw = -nw -1
        s = -1
    end
    iw = get_sign_iw(nw,Ngamma)
    return s*gamma[x,iw]
end

function get_sign_iw(nw::Integer,N::Integer)
    s = sign(nw)
    nw_bounds = min(abs(nw), N-1)
    return s*nw_bounds+1
end


# 2 particle Vertices
@inline function convertFreqArgs(ns,nt,nu,Nw)
    # assert ns+nt+nu odd, uses pos freq <-> neg freq symmetry, cuts off each freq at (Nw-1)/3
    swaps1 = ns*nu < 0
    swaps2 = nt*nu < 0
    ns,nt,nu = abs.((ns,nt,nu))
    ns = min( ns, Nw - 1 - (ns+Nw-1)%2)
    nt = min( nt, Nw - 1 - (nt+Nw-1)%2)
    nu = min( nu, Nw - 1 - (nu+Nw-1)%2)

    return ns,nt,nu,swaps1,swaps2
end

# flavor 1:xxxx, 2:zzzz, 3:xxyy, 4:xxzz, 5:zzxx, 6:xyxy, 7:xzxz, 8:zxzx
# Even though we could define less vertices, we define 8, in order to be able to use symmetries (Majorana exchange & time reversal), to only compute positive Matsubara frequencies
# Some Majorana exchanges change the flavor of the vertex (e.g. -t -> t is achieved by 1234 -> 3412 [+ time reversal] and that changes xxzz -> zzxx)
@inline function V_(Γ::VertexTypeXXZ, flavor::Int64, Rij::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    # @assert (ns+nt+nu) %2 != 0
    ns,nt,nu,swaps1,swaps2 = convertFreqArgs(ns,nt,nu,N)
    Rij = ifelse(swaps2,Rji,Rij)
    if swaps2 && (flavor == 4 || flavor == 5)
        flavor = 9-flavor
    elseif swaps1 && (flavor == 7 || flavor == 8)
        flavor = 15-flavor
    end

    return getfield(Γ,flavor)[Rij,ns+1,nt+1,nu+1]
end