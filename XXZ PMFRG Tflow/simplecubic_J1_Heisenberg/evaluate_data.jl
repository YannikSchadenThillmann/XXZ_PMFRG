
## plot structure factor
using HDF5, SpinFRGLattices, CairoMakie, FRGLatticeEvaluation, StaticArrays, Dierckx
include("/Users/yannikschaden/Julia/Helpers/misc.jl")

file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/cluster/PMFRG_XXZ_Tflow_Cubic_NLen=10_J1=1.0_N=48_NLen=10_theta=0.25pi.h5") # cluster AFM
#file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/local/PMFRG_XXZ_Tflow_Cubic_NLen=6_J1=1.414_N=8_theta=0.25pi.h5") # local AFM
#file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/local/localPMFRG_XXZ_Tflow_Cubic_NLen=6_J1=-1.4142135623730951_N=8_NLen=6_theta=0.25pi.h5") # local FM
NLen = read(file,"Par/System/NLen")
System = SimpleCubic.getCubic(NLen)
Lattice = LatticeInfo(System,SimpleCubic)
Nsites = length(Lattice.SiteList)
println("Nsites = ",Nsites)
Npairs = read(file,"Par/System/Npairs")
Tvals = read(file,"State/Tvals")
Tposis = [30,80,131]

for Tpos in Tposis
    println(Tvals[Tpos])

    chi = read(file,"State/chizz")
    chi = [chi[i,1,j] + 2*sum(chi[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)] # equal tau

    chi_F = getFourier(chi[:,Tpos], Lattice)
    Q = [pi,pi,pi]
    qa = [1,0,0]

    # kspace and ktilde
    Nk = 501
    Nk2 = Int(ceil(Nk/2))
    Nk4 = Int(ceil(Nk2/2))
    kvals = LinRange(-2pi,2pi,Nk)
    ktvals = kvals/pi # ktilde values

    # BZ
    b1 = [1,0]
    b2 = [0,1]
    BZ = [b1+b2,b1-b2,-b1-b2,-b1+b2,b1+b2]
    xBZ = [el[1] for el in BZ]
    yBZ = [el[2] for el in BZ]

    chi_k = [chi_F(k1,k2,0.0) for k1 in kvals, k2 in kvals]
    min = Min(chi_k)
    chi_k = chi_k .- min
    max = Max(chi_k)
    println("max = ",max)

    fig = Figure(fontsize = 32)
    ax = Axis(fig[1,1], xlabel = L"(h,0,0)", ylabel = L"(0,l,0)", xticks = ([-1,0,1],["-1","0","1"]), yticks = ([-1,0,1],["-1","0","1"]), aspect = 1)
    hm = heatmap!(ax,ktvals,ktvals,chi_k/max, colorrange = [0.0,1.0])
    Colorbar(fig[1,2],hm, halign = :left, height = Relative(1.0))

    lines!(ax,xBZ,yBZ, color = :white, linewidth = 0.5)

    display("image/png",fig)
end

## path through BZ
fig = Figure(fontsize = 32)
ax = Axis(fig[1,1], xlabel = L"(h,h,0)", ylabel = L"(0,l,0)", xticks = ([-1,0,1],["-1","0","1"]), yticks = ([-1,0,1],["-1","0","1"]), aspect = 1)

Tposis = [30,80,131]

max = 0.0
for Tpos in Tposis
    chi = read(file,"State/chizz")
    chi = [chi[i,1,j] + 2*sum(chi[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)] # equal tau

    chi_F = getFourier(chi[:,Tpos], Lattice)
    chi_k = [chi_F(k1,k1,0.0) for k1 in kvals]
    min = Min(chi_k)
    chi_0 = chi_F(1pi,1pi,0.0) - min
    if chi_0 > max
        max = chi_0
    end
end

for Tpos in Tposis
    chi = read(file,"State/chizz")
    chi = [chi[i,1,j] + 2*sum(chi[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)] # equal tau

    chi_F = getFourier(chi[:,Tpos], Lattice)
    chi_k = [chi_F(k1,k1,0.0) for k1 in kvals]
    min = Min(chi_k)
    chi_k = chi_k .- min
    println(Max(chi_k)/max)

    lines!(ax,ktvals,chi_k/max, label = L"T = %$(strno(Tvals[Tpos],2))J^z")
end

axislegend(ax)
display("image/png",fig)








## determine Tc values from correlation length (Cluster data)
using HDF5, SpinFRGLattices, CairoMakie, FRGLatticeEvaluation, StaticArrays, Dierckx

folder = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/cluster"
filenames = readdir(folder)
filter!(x->x!=".DS_Store",filenames)
filter!(x->x!="data_analysis",filenames)

NLenthetavals = []
for filename in filenames
    file = h5open(folder*"/"*filename)
    theta = read(file,"Par/NumericalParams/theta")
    NLen = read(file,"Par/System/NLen")
    push!(NLenthetavals,(NLen,theta))
end

NLenthetavals = [round.(el, digits = 9) for el in NLenthetavals]
thetavals = unique!([el[2] for el in NLenthetavals])
NLenvals = Int.(unique!([el[1] for el in NLenthetavals]))

xi_T_vals_arr = []
Tvals_arr = []
Qqa_arr = []
iextra = []
Tcextra = Float64[]
Tcvalsp = Vector{Vector{Float64}}()

CorrelationLength2(Chi, Q::AbstractVector, qa::AbstractVector) = 1 / norm(Q - qa) * sqrt(abs(Chi(Q) / Chi(qa) - 1))
function CorrelationLength2(Chi::AbstractArray,Q::AbstractVector,direction::AbstractVector,Lattice::AbstractLattice) 
    ChiFunc = getFourier(Chi,Lattice)
    Nsites = length(Lattice.SiteList)
    L = (3*Nsites/(4pi))^(1/3)
    global qa = Q+direction/norm(direction) * 2pi/Lattice.System.NLen
    CorrelationLength2(ChiFunc,Q,qa)
end

fig = Figure()
ax = Axis(fig[1,1], title = L"J^x = sin(θ), J^z = cos(θ)", xlabel = L"\text{test}")

xlims!(ax, -0.6,1.6)
ylims!(ax, -0.05,1.1)
scatter!(ax, thetavals/pi, zeros(length(thetavals)), markersize = 8)

# determine Tc values
sphere(r,θ,ϕ) = [r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ)]

# compute xi(T)/NLen for different NLen
Tvals_f = 0.005:0.001:1.5

Nsites_vals = []
for theta in thetavals
    xi_T_vals = []
    L_vals = []
    Tvals_loc = []
    Qqa = []

    for NLen in NLenvals
        file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/cluster/PMFRG_XXZ_Tflow_Cubic_NLen=$(NLen)_J1=1.0_N=48_NLen=$(NLen)_theta=0.25pi.h5")
        key = keys(file)
        state = filter!(x->x≠"Par", key)[1]
        Tvals = read(file,state*"/Tvals")
        Tvals = filter(x->x>0.005, Tvals)
        theta = read(file,"Par/NumericalParams/theta")

        # look for Q (maximum of structure factor) at T = 0.1J or smallest value
        Tpos = MinPos(abs.(Tvals.-0.1))[1]
        
        Npairs = read(file,"Par/System/Npairs")
        
        state = read(file,state)
        chixx_R_nu_T = state["chixx"]
        chizz_R_nu_T = state["chizz"]
        
        # compute equal tau structure factor
        chixx_R_T = [chixx_R_nu_T[i,1,j] + 2*sum(chixx_R_nu_T[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)]
        chizz_R_T = [chizz_R_nu_T[i,1,j] + 2*sum(chizz_R_nu_T[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)]
        
        System = SimpleCubic.getCubic(NLen)
        Lattice = LatticeInfo(System,SimpleCubic)
        
        Nsites = length(Lattice.SiteList)
        L = (3*Nsites/(4pi))^(1/3)
        println(Nsites)
    
        # compute correlation length
        chixx_Fourier = getFourier(chixx_R_T[:,Tpos], Lattice)
        Q = [pi,pi,pi]
        qa = [1,0,0]
        xixx_T = [CorrelationLength2(chixx_R_T[:,i],Q,qa,Lattice) for i in eachindex(Tvals)]
        push!(Qqa,[Q,qa])

        chizz_Fourier = getFourier(chizz_R_T[:,Tpos], Lattice)
        Q = [pi,pi,pi]
        qa = [1,0,0]
        xizz_T = [CorrelationLength2(chizz_R_T[:,i],Q,qa,Lattice) for i in eachindex(Tvals)]
        push!(Qqa,[Q,qa])
        
        push!(xi_T_vals, (xixx_T/L,xizz_T/L))
        push!(Nsites_vals, Nsites)
        push!(Tvals_loc,round.(Tvals, digits = 6))
    end
    
    push!(xi_T_vals_arr, xi_T_vals)
    push!(Tvals_arr, Tvals_loc)
    push!(Qqa_arr,Qqa)
end

# Tc values from raw xi_T/L, xi_T_vals_arr[itheta][iNLen][xx/zz]
function xcrossingsLinInterp(arr1,arr2,xvals,points)
    dT = xvals[points].-xvals[points.+1]
    dy1 = arr2[points.+1].-arr1[points.+1]
    dy2 = arr2[points].-arr1[points]
    crossings = xvals[points.+1] .+ dT.*dy1./(dy1.-dy2)

    return crossings
end

function ycrossingsLinInterp(arr1,arr2,points)
    dy1 = arr2[points.+1].-arr1[points.+1]
    dy2 = arr2[points].-arr1[points]
    dy3 = arr2[points].-arr2[points.+1] 
    crossings = arr2[points.+1] .+ dy3.*dy1./(dy1.-dy2)

    return crossings
end

test2 = Vector{Float64}()
for mu in 1:2
    for itheta in eachindex(thetavals)
        # filter Tvals, which are available for all NLen
        Tvals = Tvals_arr[itheta][1]
        for i in 1:3
            Tvals = filter(x->x in Tvals_arr[itheta][i], Tvals)
        end
    
        pos1 = findall(x->x in Tvals, Tvals_arr[itheta][1])
        pos2 = findall(x->x in Tvals, Tvals_arr[itheta][2])
        pos3 = findall(x->x in Tvals, Tvals_arr[itheta][3])
    
        # detect crossings of xi(T)
        xi_T_vals_1 = xi_T_vals_arr[itheta][1][mu][pos2]
        xi_T_vals_2 = xi_T_vals_arr[itheta][2][mu][pos2]
        xi_T_vals_3 = xi_T_vals_arr[itheta][3][mu][pos3]
        signfunc12 = (xi_T_vals_1[1:end-1]-xi_T_vals_2[1:end-1]) .* (xi_T_vals_1[2:end]-xi_T_vals_2[2:end])
        signchanges12 = findall(x->x<0,signfunc12)
        signfunc23 = (xi_T_vals_3[1:end-1]-xi_T_vals_2[1:end-1]) .* (xi_T_vals_3[2:end]-xi_T_vals_2[2:end])
        signchanges23 = findall(x->x<0,signfunc23)
        signfunc31 = (xi_T_vals_3[1:end-1]-xi_T_vals_1[1:end-1]) .* (xi_T_vals_3[2:end]-xi_T_vals_1[2:end])
        signchanges31 = findall(x->x<0,signfunc31)
    
        # linear interpolation between two points to get Tc
        Tcs12 = xcrossingsLinInterp(xi_T_vals_1,xi_T_vals_2,Tvals,signchanges12)
        Tcs23 = xcrossingsLinInterp(xi_T_vals_2,xi_T_vals_3,Tvals,signchanges23)
        Tcs31 = xcrossingsLinInterp(xi_T_vals_3,xi_T_vals_1,Tvals,signchanges31)

        # check if collapse happens
        maxxi = maximum(union(xi_T_vals_1,xi_T_vals_2,xi_T_vals_3))
        Ttolerance = 0.1
        xitolerance = 0.2
        Tcs = Float64[]

        for i in eachindex(Tcs12)
            for j in eachindex(Tcs23)
                if abs(Tcs12[i]-Tcs23[j]) < Ttolerance
                    for k in eachindex(Tcs31)
                        if abs(Tcs12[i]-Tcs31[k]) < Ttolerance && abs(Tcs31[k]-Tcs23[j]) < Ttolerance
                            y12 = ycrossingsLinInterp(xi_T_vals_1,xi_T_vals_2,signchanges12)
                            y23 = ycrossingsLinInterp(xi_T_vals_2,xi_T_vals_3,signchanges23)
                            y31 = ycrossingsLinInterp(xi_T_vals_3,xi_T_vals_1,signchanges31)
                
                            if itheta == 81 || abs(y12[i]-y23[j]) < xitolerance*maxxi && abs(y23[j]-y31[k]) < xitolerance*maxxi && abs(y31[k]-y12[i]) < xitolerance*maxxi
                                push!(Tcs,Tcs23[j])
                            end
                        end
                    end
                end
            end
        end 

        if length(Tcs) > 0
            push!(Tcvalsp, Tcs)
        else
            push!(Tcvalsp, Tcs)
        end
    end
end

Tcvalsp

Ntheta = length(thetavals)
Tcvalspxxp = Tcvalsp[1:Ntheta]
Tcvalspzzp = Tcvalsp[Ntheta+1:end]

Tc = Tcvalspxxp[1]*sqrt(2)

## plot xi
fig = Figure(fontsize = 16)
ax = Axis(fig[1,1], xlabel = L"T/J^z", ylabel = L"\xi/L", aspect = 1.8)
scatter!(ax,[1.0],[2.0],color = :white, label = L"N_\text{sites}")
NLenvalsp = NLenvals[1:end-1]

xi_14_15 = []
for iN in eachindex(NLenvalsp)
    xvals = sqrt(2)*Tvals_arr[1][iN]
    yvals = xi_T_vals_arr[1][iN][1]
    spline = Spline1D(reverse(xvals),reverse(yvals))
    
    scatter!(ax,xvals,yvals, label = L"%$(Nsites_vals[iN])", markersize = 8)
    xvals = sqrt(2)*Tvals_arr[1][iN][end]:0.01:2.0
    lines!(ax,xvals,spline(xvals))
    Tc = 0.9724710140865505
end
lines!(ax,[Tc,Tc],[-0.5,2.0], color = :red, label = L"T_c/J^z")

xlims!(ax,0.89,1.06)
ylims!(ax,0.26,0.51)
axislegend(ax)
display("image/png",fig)


## staggered magnetization
using HDF5, SpinFRGLattices, CairoMakie, FRGLatticeEvaluation, StaticArrays, Dierckx

file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/cluster/PMFRG_XXZ_Tflow_Cubic_NLen=14_J1=1.0_N=48_NLen=14_theta=0.25pi.h5")
NLen = read(file,"Par/System/NLen")
System = SimpleCubic.getCubic(NLen)
Lattice = LatticeInfo(System,SimpleCubic)
Nsites = length(Lattice.SiteList)
Npairs = read(file,"Par/System/Npairs")
Tvals = read(file,"State/Tvals")

chi = read(file,"State/chizz")
chi = [chi[i,1,j] + 2*sum(chi[i,2:end,j]) for i in 1:Npairs, j in eachindex(Tvals)] # equal tau

m_arr = System.siteSum.m[:,1] # multiplicity of site occuring for each inequivalent pair
deleteat!(m_arr, m_arr .== 0)

chistaggered = [sum([abs(chi[i,Tj]) * m_arr[i] for i in 1:Npairs]) for Tj in eachindex(Tvals)]/Nsites
