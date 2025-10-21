
# magnetic Susceptibility
function V2_(Vertex::AbstractArray, Rij::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns,nt,nu,swapsites = convertFreqArgs(ns,nt,nu,N)
    Rij = ifelse(swapsites,Rji,Rij)
    return Vertex[Rij,ns+1,nt+1,nu+1]
end

# susceptibility with size(χ) = (Npairs, N, T) [Matsubara frequency resolved for equal τ correlator]
function getChi1(flavor::Int64, gammax::AbstractArray, gammaz::AbstractArray, Γ::AbstractArray, T::Real, Par, Numax)
	(;N,lenIntw_acc,np_vec) = Par.NumericalParams
	(;Npairs,invpairs,PairTypes,OnsitePairs) = Par.System

	@inline iGx(x,w) = iG_(gammax,x,T,w)
    @inline iGz(x,w) = iG_(gammaz,x,T,w)
    @inline Γ_(Rij,s,t,u) = V2_(Γ,Rij,s,t,u,invpairs[Rij],N)

	Chi = zeros(_getFloatType(Par), Npairs, N)

	@inbounds Threads.@threads for Rij in 1:Npairs
		(;xi,xj) = PairTypes[Rij]
		for i_nu in 1:Numax
       		n_nu = np_vec[i_nu]
		
			for nw1 in -lenIntw_acc:lenIntw_acc-1
				if Rij in OnsitePairs
                    if flavor == 1
                        Chi[Rij,i_nu] += iGx(xi,nw1) *iGz(xi,nw1+n_nu)
                    elseif flavor == 3
                        Chi[Rij,i_nu] += iGx(xi,nw1) *iGx(xi,nw1+n_nu)
                    end
				end
				for nw2 in -lenIntw_acc:lenIntw_acc-1
					npw1pw2 = n_nu+nw1+nw2+1
					w2mw1 = nw2-nw1
                    if flavor == 1
                        GGGG = iGx(xi,nw2+n_nu) *iGz(xi,nw2) *iGx(xj,nw1) *iGz(xj,nw1+n_nu)
                    elseif flavor == 3
                        GGGG = iGx(xi,nw2+n_nu) *iGx(xi,nw2) *iGx(xj,nw1) *iGx(xj,nw1+n_nu)
                    end
                    Chi[Rij,i_nu] += GGGG *Γ_(Rij,0,npw1pw2,w2mw1)
                end
            end
        end
    end
	return(Chi)
end
getChixx(State::ArrayPartition, T::Real,Par, Numax) = getChi1(1,State.x[2],State.x[3],State.x[10],T,Par, Numax)
getChizz(State::ArrayPartition, T::Real,Par, Numax) = getChi1(3,State.x[2],State.x[3],State.x[9],T,Par, Numax)

getChixx(State::AbstractVector, T::Real,Par, Numax) = getChi1(1,State[2],State[3],State[10],T,Par, Numax)
getChizz(State::AbstractVector, T::Real,Par, Numax) = getChi1(3,State[2],State[3],State[9],T,Par, Numax)
