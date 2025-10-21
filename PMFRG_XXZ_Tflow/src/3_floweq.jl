
# define r.h.s of flow equations
function getDeriv!(Deriv,State,setup,T, Katanin)
    (;X,Par,VertexBuffers) = setup
    Workspace = OneLoopWorkspaceXXZ(Deriv,State,X,Par)

    getDFint!(Workspace,T)
    get_Self_Energy!(Workspace,T)
    getXBubble!(Workspace,T,VertexBuffers, Katanin)
    symmetrizeBubble!(Workspace.X,Par)
    addToVertexFromBubble!(Workspace.Deriv.Γ, Workspace.X)
    #symmetrizeVertex!(Workspace.Deriv.Γ, Par)
    return
end

# free energy
function getDFint!(Workspace,T::Real)
    (;State,Deriv,Par) = Workspace 
    (;lenIntw_acc) = Par.NumericalParams 
    NUnique = Par.System.NUnique 
    
    @inline γx(xk,nw) = gamma_(State.γ.x, xk, nw)
    @inline γz(xk,nw) = gamma_(State.γ.z, xk, nw)
    @inline iGx(xk,nw) = iG_(State.γ.x, xk, T, nw)
    @inline iGz(xk,nw) = iG_(State.γ.z, xk, T, nw)
    @inline iSx(xk,nw) = iS_(State.γ.x, xk, T, nw)
    @inline iSz(xk,nw) = iS_(State.γ.z, xk, T, nw)
    
    for xk in 1:NUnique
        sumres = 0.
        for nw in -lenIntw_acc:lenIntw_acc-1
            sumres +=
            (
                +iSx(xk,nw)/iGx(xk,nw) *iG0(nw,T) *γx(xk,nw) *2
                +iSz(xk,nw)/iGz(xk,nw) *iG0(nw,T) *γz(xk,nw)
            )
        end
        Deriv.f_int[xk] = -1/2*sumres
    end
end

# self energy
function get_Self_Energy!(Workspace, T)
    Par = Workspace.Par
    @inline iSx(xk,nw) = iS_(Workspace.State.γ.x, xk, T, nw)/2
    @inline iSz(xk,nw) = iS_(Workspace.State.γ.z, xk, T, nw)/2
    compute1PartBubble!(Workspace.Deriv.γ.x,Workspace.Deriv.γ.z, Workspace.State.Γ, iSx, iSz, Par)
end

function compute1PartBubble!(Dgammax::AbstractArray,Dgammaz::AbstractArray,Γ::VertexTypeXXZ,Px::Function,Pz::Function,Par)
    invpairs = Par.System.invpairs

    Dgammax = setZero!(Dgammax)
    Dgammaz = setZero!(Dgammaz)

    @inline Γxxxx(Rij,s,t,u) = V_(Γ,1,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)
    @inline Γzzzz(Rij,s,t,u) = V_(Γ,2,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)

    @inline Γxxyy(Rij,s,t,u) = V_(Γ,3,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)
    @inline Γxxzz(Rij,s,t,u) = V_(Γ,4,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)
    @inline Γzzxx(Rij,s,t,u) = V_(Γ,5,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)

    addTo1PartBubble!(Dgammax,Dgammaz,Γxxxx,Γzzzz,Γxxyy,Γxxzz,Γzzxx,Px,Pz,Par)
end

function addTo1PartBubble!(Dgammax::AbstractArray,Dgammaz::AbstractArray,Γxxxx::Function,Γzzzz::Function,Γxxyy::Function,Γxxzz::Function,Γzzxx::Function,Px::Function,Pz::Function,Par)

    (;Ngamma,lenIntw_acc,np_vec_gamma) = Par.NumericalParams
    (;siteSum,Nsum,OnsitePairs) = Par.System

    Threads.@threads for iw1 in 1:Ngamma
        nw1 = np_vec_gamma[iw1]
        for (x,Rx) in enumerate(OnsitePairs)
            for nw in -lenIntw_acc:lenIntw_acc-1
                jsumx = 0.
                jsumz = 0.
                wpw1 = nw1+nw+1
                wmw1 = nw-nw1
                for k_spl in 1:Nsum[Rx]
                    (;m,ki,xk) = siteSum[k_spl,Rx]
                    jsumx += (Px(xk,nw)*(Γxxxx(ki,0,wmw1,wpw1) + Γxxyy(ki,0,wmw1,wpw1)) + Pz(xk,nw)*Γzzxx(ki,0,wmw1,wpw1))*m # flow equations for self energy
                    jsumz += (2*Px(xk,nw)*Γxxzz(ki,0,wmw1,wpw1) + Pz(xk,nw)*Γzzzz(ki,0,wmw1,wpw1))*m
                end
                Dgammax[x,iw1] += -jsumx
                Dgammaz[x,iw1] += -jsumz
            end
        end
    end
    return Dgammax,Dgammaz
end


# vertices
function mixedFrequencies(ns,nt,nu,nwpr)
    nw1=Int((ns+nt+nu-1)/2)
    nw2=Int((ns-nt-nu-1)/2)
    nw3=Int((-ns+nt-nu-1)/2)
    nw4=Int((-ns-nt+nu-1)/2)
    wpw1 = nwpr + nw1 + 1
    mwmw2 = - nwpr - nw2 - 1
    mwpw3 = - nwpr + nw3
    mwpw4 = - nwpr + nw4
    # @assert (ns + wmw3 +wmw4)%2 != 0 "error in freq"
    return wpw1,mwmw2,mwpw3,mwpw4
end

function getXBubble!(Workspace,T,VertexBuffers, Katanin)
    Par = Workspace.Par
    (;N,lenIntw,np_vec) = Par.NumericalParams
    (;NUnique) = Par.System
    
    @inline iG_x(xk,nw) = iG_(Workspace.State.γ.x, xk, T, nw)
    @inline iG_z(xk,nw) = iG_(Workspace.State.γ.z, xk, T, nw)
    @inline Prop_x(xk,nw) = iSKat_(Workspace.State.γ.x, Workspace.Deriv.γ.x, xk, T, nw, Katanin)
    @inline Prop_z(xk,nw) = iSKat_(Workspace.State.γ.z, Workspace.Deriv.γ.z, xk, T, nw, Katanin)

    function getKataninPropX!(BubbleProp,nw1,nw2)
        for i in 1:Par.System.NUnique, j in 1:Par.System.NUnique
            BubbleProp[1,i,j] = Prop_x(i,nw1) *iG_x(j,nw2)
            BubbleProp[2,i,j] = Prop_z(i,nw1) *iG_z(j,nw2)
            BubbleProp[3,i,j] = Prop_x(i,nw1) *iG_z(j,nw2)
            BubbleProp[4,i,j] = Prop_z(i,nw1) *iG_x(j,nw2)
        end
        return BubbleProp
    end

    function getKataninPropXTilde!(BubbleProp,nw1,nw2)
        for i in 1:Par.System.NUnique, j in 1:Par.System.NUnique
            BubbleProp[1,i,j] = (Prop_x(i,nw1) *iG_x(j,nw2) + iG_x(i,nw1) *Prop_x(j,nw2))
            BubbleProp[2,i,j] = (Prop_z(i,nw1) *iG_z(j,nw2) + iG_z(i,nw1) *Prop_z(j,nw2))
            BubbleProp[3,i,j] = (Prop_x(i,nw1) *iG_z(j,nw2) + iG_x(i,nw1) *Prop_z(j,nw2))
            BubbleProp[4,i,j] = (Prop_z(i,nw1) *iG_x(j,nw2) + iG_z(i,nw1) *Prop_x(j,nw2))
        end
        return BubbleProp
    end

    @sync for is in 1:N,it in 1:N
        Threads.@spawn begin
            Buffer = take!(VertexBuffers)

            BubblePropX = zeros(4,NUnique,NUnique)
            BubblePropXTilde = zeros(4,NUnique,NUnique)
            ns = np_vec[is]
            nt = np_vec[it]
            for nw in -lenIntw:lenIntw-1 # Matsubara sum
                MatrixPropsX = getKataninPropX!(BubblePropX,nw,nw+ns)
                MatrixPropsXTilde = getKataninPropXTilde!(BubblePropXTilde,nw,nw+ns)
                for iu in 1:N
                    nu = np_vec[iu]
                    if (ns+nt+nu)%2 == 0 # skip unphysical bosonic frequency combinations
                        continue
                    end
                    addXTilde!(Workspace,is,it,iu,nw,MatrixPropsXTilde) # add to XTilde-type bubble functions
                    if(!Par.Options.usesymmetry || nu<=nt)
                        addX!(Workspace,is,it,iu,nw,MatrixPropsX,Buffer) # add to X-type bubble functions
                    end
                end
            end
            put!(VertexBuffers, Buffer)

        end
    end
end

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, MatrixProps,VBuffer::VertexBufferType)
    (;State,X,Par) = Workspace 
    (;N,np_vec) = Par.NumericalParams
    (;Npairs,Nsum,siteSum,invpairs) = Par.System

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, mwmw2, mwpw3, mwpw4 = mixedFrequencies(ns,nt,nu,nwpr)

    bufferV_!(VBuffer.V_xxxx_21, State.Γ,1, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_zzzz_21, State.Γ,2, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_xxyy_21, State.Γ,3, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_xxzz_21, State.Γ,4, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_zzxx_21, State.Γ,5, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_xyxy_21, State.Γ,6, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_xzxz_21, State.Γ,7, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_zxzx_21, State.Γ,8, ns,mwmw2,wpw1,invpairs,N)
    bufferV_!(VBuffer.V_xyxy_12, State.Γ,6, ns,wpw1,mwmw2,invpairs,N)
    bufferV_!(VBuffer.V_xzxz_12, State.Γ,7, ns,wpw1,mwmw2,invpairs,N)
    bufferV_!(VBuffer.V_zxzx_12, State.Γ,8, ns,wpw1,mwmw2,invpairs,N)
    bufferV_!(VBuffer.V_xxxx_34, State.Γ,1, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_zzzz_34, State.Γ,2, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_xxyy_34, State.Γ,3, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_xxzz_34, State.Γ,4, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_zzxx_34, State.Γ,5, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_xyxy_34, State.Γ,6, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_xzxz_34, State.Γ,7, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_zxzx_34, State.Γ,8, ns,mwpw3,mwpw4,invpairs,N)
    bufferV_!(VBuffer.V_xyxy_43, State.Γ,6, ns,mwpw4,mwpw3,invpairs,N)
    bufferV_!(VBuffer.V_xzxz_43, State.Γ,7, ns,mwpw4,mwpw3,invpairs,N)
    bufferV_!(VBuffer.V_zxzx_43, State.Γ,8, ns,mwpw4,mwpw3,invpairs,N)

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    @inbounds for Rij in 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        X_xxxx_sum = 0.
        X_zzzz_sum = 0.

        X_xxyy_sum = 0.
        X_xxzz_sum = 0.
        X_zzxx_sum = 0.

        X_xyxy_sum = 0.
        X_xzxz_sum = 0.
        X_zxzx_sum = 0.
        @inbounds for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ik, kj, m, xk = S_ki[k_spl,Rij], S_kj[k_spl,Rij], S_m[k_spl,Rij], S_xk[k_spl,Rij]
            Pxx = MatrixProps[1,xk,xk]*m
            Pzz = MatrixProps[2,xk,xk]*m
            Pxz = MatrixProps[3,xk,xk]*m
            Pzx = MatrixProps[4,xk,xk]*m

            # prepare vertices
            # 21 structure (ik sites)
            V_xxxx_21 = VBuffer.V_xxxx_21[ik]
            V_zzzz_21 = VBuffer.V_zzzz_21[ik]

            V_xxyy_21 = VBuffer.V_xxyy_21[ik]
            V_xxzz_21 = VBuffer.V_xxzz_21[ik]
            V_zzxx_21 = VBuffer.V_zzxx_21[ik]

            V_xyxy_21 = VBuffer.V_xyxy_21[ik]
            V_xzxz_21 = VBuffer.V_xzxz_21[ik]
            V_zxzx_21 = VBuffer.V_zxzx_21[ik]

            # 12 structure (ik sites)
            V_xyxy_12 = VBuffer.V_xyxy_12[ik]
            V_xzxz_12 = VBuffer.V_xzxz_12[ik]
            V_zxzx_12 = VBuffer.V_zxzx_12[ik]

            # 34 structure (kj sites)
            V_xxxx_34 = VBuffer.V_xxxx_34[kj]
            V_zzzz_34 = VBuffer.V_zzzz_34[kj]

            V_xxyy_34 = VBuffer.V_xxyy_34[kj]
            V_xxzz_34 = VBuffer.V_xxzz_34[kj]
            V_zzxx_34 = VBuffer.V_zzxx_34[kj]

            V_xyxy_34 = VBuffer.V_xyxy_34[kj]
            V_xzxz_34 = VBuffer.V_xzxz_34[kj]
            V_zxzx_34 = VBuffer.V_zxzx_34[kj]

            # 32 structure (kj sites)
            V_xyxy_43 = VBuffer.V_xyxy_43[kj]
            V_xzxz_43 = VBuffer.V_xzxz_43[kj]
            V_zxzx_43 = VBuffer.V_zxzx_43[kj]

            # flow equations
            X_xxxx_sum += -Pxx *(V_xxxx_21 *V_xxxx_34 + V_xxyy_21 *V_xxyy_34) - Pzz *(V_xxzz_21 *V_zzxx_34)
            X_zzzz_sum += -Pxx *(V_zzxx_21 *V_xxzz_34 + V_zzxx_21 *V_xxzz_34) - Pzz *(V_zzzz_21 *V_zzzz_34) # first two terms equal

            X_xxyy_sum += -Pxx *(V_xxxx_21 *V_xxyy_34 + V_xxyy_21 *V_xxxx_34) - Pzz *(V_xxzz_21 *V_zzxx_34)
            X_xxzz_sum += -Pxx *(V_xxxx_21 *V_xxzz_34 + V_xxyy_21 *V_xxzz_34) - Pzz *(V_xxzz_21 *V_zzzz_34)
            X_zzxx_sum += -Pxx *(V_zzxx_21 *V_xxxx_34 + V_zzxx_21 *V_xxyy_34) - Pzz *(V_zzzz_21 *V_zzxx_34)

            X_xyxy_sum +=  Pxx *V_xyxy_12 *V_xyxy_34 + Pxx *V_xyxy_21 *V_xyxy_43
            X_xzxz_sum +=  Pxz *V_xzxz_12 *V_xzxz_34 + Pzx *V_xzxz_21 *V_zxzx_43
            X_zxzx_sum +=  Pxz *V_zxzx_21 *V_xzxz_43 + Pzx *V_zxzx_12 *V_zxzx_34
        end
        X.xxxx[Rij,is,it,iu] += X_xxxx_sum
        X.zzzz[Rij,is,it,iu] += X_zzzz_sum

        X.xxyy[Rij,is,it,iu] += X_xxyy_sum
        X.xxzz[Rij,is,it,iu] += X_xxzz_sum
        X.zzxx[Rij,is,it,iu] += X_zzxx_sum

        X.xyxy[Rij,is,it,iu] += X_xyxy_sum
        X.xzxz[Rij,is,it,iu] += X_xzxz_sum
        X.zxzx[Rij,is,it,iu] += X_zxzx_sum
    end
    return
end

function addXTilde!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, MatrixProps)
    (;State,X,Par) = Workspace 
    (;N,np_vec) = Par.NumericalParams
    (;Npairs,invpairs,PairTypes,OnsitePairs) = Par.System

    @inline Vxxxx_(Rij,s,t,u) = V_(State.Γ,1,Rij,s,t,u,invpairs[Rij],N)
    @inline Vzzzz_(Rij,s,t,u) = V_(State.Γ,2,Rij,s,t,u,invpairs[Rij],N)

    @inline Vxxyy_(Rij,s,t,u) = V_(State.Γ,3,Rij,s,t,u,invpairs[Rij],N)
    @inline Vxxzz_(Rij,s,t,u) = V_(State.Γ,4,Rij,s,t,u,invpairs[Rij],N)
    @inline Vzzxx_(Rij,s,t,u) = V_(State.Γ,5,Rij,s,t,u,invpairs[Rij],N)

    @inline Vxyxy_(Rij,s,t,u) = V_(State.Γ,6,Rij,s,t,u,invpairs[Rij],N)
    @inline Vxzxz_(Rij,s,t,u) = V_(State.Γ,7,Rij,s,t,u,invpairs[Rij],N)
    @inline Vzxzx_(Rij,s,t,u) = V_(State.Γ,8,Rij,s,t,u,invpairs[Rij],N)

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, mwmw2, mwpw3, mwpw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #the following Xtilde need to be computed (and can be computed) only for nonlocal pairs Rij != Rii, because XTilde^ii = X^ii
    for Rij in 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        (;xi,xj) = PairTypes[Rij]
        # @views (Pxx,Pzz,Pxz,Pzx) = MatrixProps[:,xi,xj] # Pxz ~ SKat_x(w) G_z(w+s) + G_x(w) SKat_z(w+s)
        Pxx = MatrixProps[1,xi,xj]
        Pzz = MatrixProps[2,xi,xj]
        Pxz = MatrixProps[3,xi,xj]
        Pzx = MatrixProps[4,xi,xj]

        # prepare vertices
        # (w1+w, s, -w-w2) structure
        V_xxxx_1s2 = Vxxxx_(Rij,wpw1,ns,mwmw2)
        V_zzzz_1s2 = Vzzzz_(Rij,wpw1,ns,mwmw2)

        V_xxyy_1s2 = Vxxyy_(Rij,wpw1,ns,mwmw2)
        V_xxzz_1s2 = Vxxzz_(Rij,wpw1,ns,mwmw2)
        V_zzxx_1s2 = Vzzxx_(Rij,wpw1,ns,mwmw2)

        V_xyxy_1s2 = Vxyxy_(Rij,wpw1,ns,mwmw2)
        V_xzxz_1s2 = Vxzxz_(Rij,wpw1,ns,mwmw2)
        V_zxzx_1s2 = Vzxzx_(Rij,wpw1,ns,mwmw2)

        # (w1+w, -w-w2, s) structure
        V_xyxy_12s = Vxyxy_(Rij,wpw1,mwmw2,ns)
        V_xzxz_12s = Vxzxz_(Rij,wpw1,mwmw2,ns)
        V_zxzx_12s = Vzxzx_(Rij,wpw1,mwmw2,ns)

        # (-w+w3, s, -w+w4) structure
        V_xxxx_3s4 = Vxxxx_(Rij,mwpw3,ns,mwpw4)
        V_zzzz_3s4 = Vzzzz_(Rij,mwpw3,ns,mwpw4)

        V_xxyy_3s4 = Vxxyy_(Rij,mwpw3,ns,mwpw4)
        V_xxzz_3s4 = Vxxzz_(Rij,mwpw3,ns,mwpw4)
        V_zzxx_3s4 = Vzzxx_(Rij,mwpw3,ns,mwpw4)

        V_xyxy_3s4 = Vxyxy_(Rij,mwpw3,ns,mwpw4)
        V_xzxz_3s4 = Vxzxz_(Rij,mwpw3,ns,mwpw4)
        V_zxzx_3s4 = Vzxzx_(Rij,mwpw3,ns,mwpw4)

        # (-w+w3, -w+w4, s) structure
        V_xyxy_34s = Vxyxy_(Rij,mwpw3,mwpw4,ns)
        V_xzxz_34s = Vxzxz_(Rij,mwpw3,mwpw4,ns)
        V_zxzx_34s = Vzxzx_(Rij,mwpw3,mwpw4,ns)

        # flow equations
        X.Txxxx[Rij,is,it,iu] += Pxx *(V_xxxx_1s2 *V_xxxx_3s4 + V_xyxy_1s2 *V_xyxy_3s4) + Pzz *V_xzxz_1s2 *V_zxzx_3s4
        X.Tzzzz[Rij,is,it,iu] += Pxx *(V_zxzx_1s2 *V_xzxz_3s4 + V_zxzx_1s2 *V_xzxz_3s4) + Pzz *V_zzzz_1s2 *V_zzzz_3s4 

        X.Txxyy[Rij,is,it,iu] += Pxx *(V_xxxx_1s2 *V_xyxy_3s4 + V_xyxy_1s2 *V_xxxx_3s4) + Pzz *V_xzxz_1s2 *V_zxzx_3s4
        X.Txxzz[Rij,is,it,iu] += Pxx *(V_xxxx_1s2 *V_xzxz_3s4 + V_xyxy_1s2 *V_xzxz_3s4) + Pzz *V_xzxz_1s2 *V_zzzz_3s4
        X.Tzzxx[Rij,is,it,iu] += Pxx *(V_zxzx_1s2 *V_xxxx_3s4 + V_zxzx_1s2 *V_xyxy_3s4) + Pzz *V_zzzz_1s2 *V_zxzx_3s4

        X.Txyxy[Rij,is,it,iu] += Pxx *V_xxyy_1s2 *V_xxyy_3s4 + Pxx *V_xyxy_12s *V_xyxy_34s
        X.Txzxz[Rij,is,it,iu] += Pxz *V_xxzz_1s2 *V_xxzz_3s4 + Pzx *V_xzxz_12s *V_zxzx_34s
        X.Tzxzx[Rij,is,it,iu] += Pzx *V_zzxx_1s2 *V_zzxx_3s4 + Pxz *V_zxzx_12s *V_xzxz_34s

        X.Txyyx[Rij,is,it,iu] += -Pxx *V_xxyy_1s2 *V_xyxy_34s - Pxx *V_xyxy_12s *V_xxyy_3s4
        X.Txzzx[Rij,is,it,iu] += -Pxz *V_xxzz_1s2 *V_xzxz_34s - Pzx *V_xzxz_12s *V_zzxx_3s4
        X.Tzxxz[Rij,is,it,iu] += -Pzx *V_zzxx_1s2 *V_zxzx_34s - Pxz *V_zxzx_12s *V_xxzz_3s4
    end
end

function symmetrizeBubble!(X::BubbleTypeXXZ,Par)
    N = Par.NumericalParams.N
    (;Npairs,OnsitePairs) = Par.System
    usesymmetry = Par.Options.usesymmetry
    # can't use the u <--> t symmetry for XXZ model
    #= if(usesymmetry)
        for it in 1:N
            for iu in it+1:N, is in 1:N, Rij in 1:Npairs
                X.a[Rij,is,it,iu] = -X.a[Rij,is,iu,it]
                X.b[Rij,is,it,iu] = -X.b[Rij,is,iu,it]
                X.c[Rij,is,it,iu] = (
                + X.a[Rij,is,it,iu]+
                - X.b[Rij,is,it,iu]+
                + X.c[Rij,is,iu,it])
            end
        end
    end =#

    #local definitions of X.Tilde vertices (needed, because local XTilde not defined)
    Threads.@threads for iu in 1:N
        for it in 1:N, is in 1:N, R in OnsitePairs
            X.Txxxx[R,is,it,iu] = X.xxxx[R,is,it,iu]
            X.Tzzzz[R,is,it,iu] = X.zzzz[R,is,it,iu]

            X.Txxyy[R,is,it,iu] = X.xxyy[R,is,it,iu]
            X.Txxzz[R,is,it,iu] = X.xxzz[R,is,it,iu]
            X.Tzzxx[R,is,it,iu] = X.zzxx[R,is,it,iu]

            X.Txyxy[R,is,it,iu] = X.xyxy[R,is,it,iu]
            X.Txzxz[R,is,it,iu] = X.xzxz[R,is,it,iu]
            X.Tzxzx[R,is,it,iu] = X.zxzx[R,is,it,iu]

            X.Txyyx[R,is,it,iu] = - X.xyxy[R,is,iu,it]
            X.Txzzx[R,is,it,iu] = - X.xzxz[R,is,iu,it]
            X.Tzxxz[R,is,it,iu] = - X.zxzx[R,is,iu,it]
        end
    end
end


# flow equations in terms of X
function addToVertexFromBubble!(Γ::VertexTypeXXZ,X::BubbleTypeXXZ)
    Threads.@threads for iu in axes(Γ.xxxx,4)
        for it in axes(Γ.xxxx,3), is in axes(Γ.xxxx,2), Rij in axes(Γ.xxxx,1)
            Γ.xxxx[Rij,is,it,iu] += X.xxxx[Rij,is,it,iu] - X.Txxxx[Rij,it,is,iu] + X.Txxxx[Rij,iu,is,it]
            Γ.zzzz[Rij,is,it,iu] += X.zzzz[Rij,is,it,iu] - X.Tzzzz[Rij,it,is,iu] + X.Tzzzz[Rij,iu,is,it]
            Γ.xxyy[Rij,is,it,iu] += X.xxyy[Rij,is,it,iu] - X.Txyxy[Rij,it,is,iu] + X.Txyxy[Rij,iu,is,it]
            Γ.xxzz[Rij,is,it,iu] += X.xxzz[Rij,is,it,iu] - X.Txzxz[Rij,it,is,iu] + X.Txzxz[Rij,iu,is,it]
            Γ.zzxx[Rij,is,it,iu] += X.zzxx[Rij,is,it,iu] - X.Tzxzx[Rij,it,is,iu] + X.Tzxzx[Rij,iu,is,it]
            Γ.xyxy[Rij,is,it,iu] += X.xyxy[Rij,is,it,iu] - X.Txxyy[Rij,it,is,iu] + X.Txyyx[Rij,iu,is,it]
            Γ.xzxz[Rij,is,it,iu] += X.xzxz[Rij,is,it,iu] - X.Txxzz[Rij,it,is,iu] + X.Txzzx[Rij,iu,is,it]
            Γ.zxzx[Rij,is,it,iu] += X.zxzx[Rij,is,it,iu] - X.Tzzxx[Rij,it,is,iu] + X.Tzxxz[Rij,iu,is,it]
        end
    end
    return Γ
end

# used in code for isotropic Heisenberg models, but not (yet) here
#= function symmetrizeVertex!(Γ::VertexTypeXXZ,Par)
    N = Par.NumericalParams.N
    for iu in 1:N
        for it in 1:N, is in 1:N, R in Par.System.OnsitePairs
            Γ.c[R,is,it,iu] = -Γ.b[R,it,is,iu]
        end
    end
end =#

@inline function bufferV_!(
    Cache,
    Γ::VertexTypeXXZ,
    flavor::Int64,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    invpairs::AbstractArray,
    N::Integer
)

    ns,nt,nu,swaps1,swaps2 = convertFreqArgs(ns,nt,nu,N)
    if swaps2 && (flavor == 4 || flavor == 5)
        flavor = 9-flavor
    elseif swaps1 && (flavor == 7 || flavor == 8)
        flavor = 15-flavor
    end

    is, it, iu = ns + 1, nt + 1, nu + 1
    
    Vertex = getfield(Γ,flavor)
    
    @inbounds begin
        if swaps2
            for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[invpairs[R], is, it, iu]
            end
        else
            for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[R, is, it, iu]
            end
        end
    end
end