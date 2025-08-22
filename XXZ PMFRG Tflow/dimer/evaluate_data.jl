
using HDF5, SpinFRGLattices, CairoMakie

folder = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/"
filenames = readdir(folder)
filenames = filter(x->x!=".DS_Store",filenames)
filenames = filter(x->x[1]!='t',filenames)

vop = false

if vop
    filenames = filter(x->x[1]=='v',filenames)
else
    filenames = filter(x->x[1]!='v',filenames)
end

Paras = []
for filename in filenames
    file = h5open(folder*filename)
    Tvals = read(file,"State/Tvals")
    theta = read(file,"Par/NumericalParams/theta")
    N = read(file,"Par/NumericalParams/N")
    NLen = read(file,"Par/System/NLen")
    Npairs = read(file,"Par/System/Npairs")
    couplings = read(file,"Par/System/couplings")

    push!(Paras,[theta,N,NLen,Npairs,couplings,Tvals])
end

thetavals = unique!([el[1] for el in Paras])
nthetavals = thetavals/pi
Npairs = unique!([el[4] for el in Paras])
Tvals = unique!([el[6] for el in Paras])

Tvals = reducer(Tvals)
Npairs = reducer(Npairs)

NT = length(Tvals)
chixx0_theta = zeros(Npairs,NT,4)
chizz0_theta = zeros(Npairs,NT,4)

for j in eachindex(filenames)
    filename = filenames[j]
    file = h5open(folder*filename)
    
    chixx0 = read(file,"State/chixx")[:,1,:]
    chizz0 = read(file,"State/chizz")[:,1,:]

    chixx0_theta[:,:,j] = chixx0
    chizz0_theta[:,:,j] = chizz0
end

fig = Figure(size = (1300,1000), title = "XXZ Dimer")
axxxloc = Axis(fig[1,1], ylabel = L"χ",xlabel = L"T",title = "χxx local", xticks = ([0.01,0.1,1.0,10.,100.],[L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"]), aspect = 4/3, xscale = log10)
axxxnonloc = Axis(fig[1,2], ylabel = L"χ",xlabel = L"T",title = "χxx non-local", xticks = ([0.01,0.1,1.0,10.,100.],[L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"]), aspect = 4/3, xscale = log10)
axzzloc = Axis(fig[2,1], ylabel = L"χ",xlabel = L"T",title = "χzz local", xticks = ([0.01,0.1,1.0,10.,100.],[L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"]), aspect = 4/3, xscale = log10)
axzznonloc = Axis(fig[2,2], ylabel = L"χ",xlabel = L"T",title = "χzz non-local", xticks = ([0.01,0.1,1.0,10.,100.],[L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"]), aspect = 4/3, xscale = log10)

models = ("AFM Ising","AFM Heisenberg","AFM XX","Jx = -Jz")

# exact susceptibilities
x(t,T) = sin(t)/(4*T)
z(t,T) = cos(t)/(4*T)
chixxloc(t,T) = ifelse(t == 0.25π, 1/(4*T) * (exp(-2*z(t,T)) + sinh(2*z(t,T))/(2*z(t,T)))/(exp(-2*z(t,T)) + cosh(2*z(t,T))), 1/(4*T) * (exp(-x(t,T)) * sinh(z(t,T)-x(t,T))/(z(t,T)-x(t,T)) + exp(x(t,T)) * sinh(z(t,T)+x(t,T))/(z(t,T)+x(t,T)))/(exp(-z(t,T)) + exp(z(t,T))*cosh(2*x(t,T))))
chixxnonloc(t,T) = ifelse(t == 0.25π, 1/(4*T) * (exp(-2*z(t,T)) - sinh(2*z(t,T))/(2*z(t,T)))/(exp(-2*z(t,T)) + cosh(2*z(t,T))), 1/(4*T) * (exp(-x(t,T)) * sinh(z(t,T)-x(t,T))/(z(t,T)-x(t,T)) - exp(x(t,T)) * sinh(z(t,T)+x(t,T))/(z(t,T)+x(t,T)))/(exp(-z(t,T)) + exp(z(t,T))*cosh(2*x(t,T))))
chizzloc(t,T) = ifelse(t == 0, 1/(4*T), 1/(4*T) * (exp(-2*z(t,T)) + sinh(2*x(t,T))/(2*x(t,T)))/(exp(-2*z(t,T))+cosh(2*x(t,T))))
chizznonloc(t,T) = ifelse(t == 0, -1/(4*T) * tanh(z(t,T)), 1/(4*T) * (exp(-2*z(t,T)) - sinh(2*x(t,T))/(2*x(t,T)))/(exp(-2*z(t,T))+cosh(2*x(t,T))))

xvals = 0.0:0.001:1.0

for j in eachindex(nthetavals)
    scatter!(axxxloc,Tvals,chixx0_theta[1,:,j], label = "XXZ code, θ = "*string(round(nthetavals[j],digits = 2))*"π ("*models[j]*")", markersize = 15)
    scatter!(axxxnonloc,Tvals,chixx0_theta[2,:,j], label = "XXZ code, θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
    scatter!(axzzloc,Tvals,chizz0_theta[1,:,j], label = "XXZ code, θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
    scatter!(axzznonloc,Tvals,chizz0_theta[2,:,j], label = "XXZ code, θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)

    t = pi * nthetavals[j]
    yvals = [chixxloc(t,x) for x in xvals]
    scatterlines!(axxxloc, xvals, yvals, linewidth = 2, markersize = 3, label = "XXZ exact, θ = "*string(round(t/pi,digits = 2))*"π")
    yvals = [chixxnonloc(t,x) for x in xvals]
    scatterlines!(axxxnonloc, xvals, yvals, linewidth = 2, markersize = 3)
    yvals = [chizzloc(t,x) for x in xvals]
    scatterlines!(axzzloc, xvals, yvals, linewidth = 2, markersize = 3)
    yvals = [chizznonloc(t,x) for x in xvals]
    scatterlines!(axzznonloc, xvals, yvals, linewidth = 2, markersize = 3)
end

# plots
xlims!(axxxloc,0.01,10.0)
xlims!(axxxnonloc,0.01,10.0)
xlims!(axzzloc,0.01,10.0)
xlims!(axzznonloc,0.01,10.0)

ylims!(axxxloc,0.0,2.0)
ylims!(axxxnonloc,-1.5,0.1)
ylims!(axzzloc,0.0,2.6)
ylims!(axzznonloc,-2.6,1.6)

#text!(axxxloc,.7,0.7; text = L"J_x = \text{sin}θ", fontsize = 24)
#text!(axxxloc,.7,0.55; text = L"J_z = \text{cos}θ", fontsize = 24)

axislegend(axxxloc)
#axislegend(axxxnonloc, position = :rb)
#axislegend(axzzloc)
#axislegend(axzznonloc, position = :rb)
fig




## linear scale x axis

using HDF5, SpinFRGLattices, CairoMakie

folder = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/"
filenames = readdir(folder)
filenames = filter(x->x!=".DS_Store",filenames)
filenames = filter(x->x[1]!='t',filenames)

vop = false

if vop
    filenames = filter(x->x[1]=='v',filenames)
else
    filenames = filter(x->x[1]!='v',filenames)
end

Paras = []
for filename in filenames
    file = h5open(folder*filename)
    Tvals = read(file,"State/Tvals")
    theta = read(file,"Par/NumericalParams/theta")
    N = read(file,"Par/NumericalParams/N")
    NLen = read(file,"Par/System/NLen")
    Npairs = read(file,"Par/System/Npairs")
    couplings = read(file,"Par/System/couplings")

    push!(Paras,[theta,N,NLen,Npairs,couplings,Tvals])
end

thetavals = unique!([el[1] for el in Paras])
nthetavals = thetavals/pi
Npairs = unique!([el[4] for el in Paras])
Tvals = unique!([el[6] for el in Paras])

Tvals = reducer(Tvals)
Npairs = reducer(Npairs)

NT = length(Tvals)
chixx0_theta = zeros(Npairs,NT,4)
chizz0_theta = zeros(Npairs,NT,4)

for j in eachindex(filenames)
    filename = filenames[j]
    file = h5open(folder*filename)
    
    chixx0 = read(file,"State/chixx")[:,1,:]
    chizz0 = read(file,"State/chizz")[:,1,:]

    chixx0_theta[:,:,j] = chixx0
    chizz0_theta[:,:,j] = chizz0
end

fig = Figure(size = (1300,1000), fontsize = 24)
axxxloc = Axis(fig[1,1], ylabel = L"χ",xlabel = L"T",title = "χxx local", aspect = 4/3)
axxxnonloc = Axis(fig[1,2], ylabel = L"χ",xlabel = L"T",title = "χxx non-local", aspect = 4/3)
axzzloc = Axis(fig[2,1], ylabel = L"χ",xlabel = L"T",title = "χzz local", aspect = 4/3)
axzznonloc = Axis(fig[2,2], ylabel = L"χ",xlabel = L"T",title = "χzz non-local", aspect = 4/3)

models = ("AFM Ising","AFM Heisenberg","AFM XX","Jx = -Jz")

# exact susceptibilities
x(t,T) = sin(t)/(4*T)
z(t,T) = cos(t)/(4*T)
chixxloc(t,T) =  1/(4*T) * (exp(-x(t,T)) * sinh(z(t,T)-x(t,T))/(z(t,T)-x(t,T)) + exp(x(t,T)) * sinh(z(t,T)+x(t,T))/(z(t,T)+x(t,T)))/(exp(-z(t,T)) + exp(z(t,T))*cosh(2*x(t,T)))
chixxnonloc(t,T) = 1/(4*T) * (exp(-x(t,T)) * sinh(z(t,T)-x(t,T))/(z(t,T)-x(t,T)) - exp(x(t,T)) * sinh(z(t,T)+x(t,T))/(z(t,T)+x(t,T)))/(exp(-z(t,T)) + exp(z(t,T))*cosh(2*x(t,T)))
chizzloc(t,T) = ifelse(t == 0, 1/(4*T), 1/(4*T) * (exp(-2*z(t,T)) + sinh(2*x(t,T))/(2*x(t,T)))/(exp(-2*z(t,T))+cosh(2*x(t,T))))
chizznonloc(t,T) = ifelse(t == 0, -1/(4*T) * tanh(z(t,T)), 1/(4*T) * (exp(-2*z(t,T)) - sinh(2*x(t,T))/(2*x(t,T)))/(exp(-2*z(t,T))+cosh(2*x(t,T))))

xvals = 0.0:0.001:1.0

for j in eachindex(nthetavals)
    scatter!(axxxloc,Tvals,chixx0_theta[1,:,j], label = "θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
    scatter!(axxxnonloc,Tvals,chixx0_theta[2,:,j], label = "θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
    scatter!(axzzloc,Tvals,chizz0_theta[1,:,j], label = "θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
    scatter!(axzznonloc,Tvals[15:end],chizz0_theta[2,:,j][15:end], label = "θ = "*string(round(nthetavals[j],digits = 2))*"π", markersize = 15)
end

for j in eachindex(nthetavals)
    t = pi * nthetavals[j]
    yvals = [chixxloc(t,x) for x in xvals]
    scatterlines!(axxxloc, xvals, yvals, linewidth = 2, markersize = 3, label = "Exact, θ = "*string(round(t/pi,digits = 2))*"π")
    yvals = [chixxnonloc(t,x) for x in xvals]
    scatterlines!(axxxnonloc, xvals, yvals, linewidth = 2, markersize = 3)
    yvals = [chizzloc(t,x) for x in xvals]
    scatterlines!(axzzloc, xvals, yvals, linewidth = 2, markersize = 3)
    yvals = [chizznonloc(t,x) for x in xvals]
    scatterlines!(axzznonloc, xvals, yvals, linewidth = 2, markersize = 3)
end

# plots
xlims!(axxxloc,0.01,0.7)
xlims!(axxxnonloc,0.01,0.7)
xlims!(axzzloc,0.01,0.7)
xlims!(axzznonloc,0.01,0.7)

ylims!(axxxloc,0.0,2.0)
ylims!(axxxnonloc,-1.5,0.1)
ylims!(axzzloc,0.0,2.6)
ylims!(axzznonloc,-2.6,1.6)

#text!(axxxloc,.7,0.7; text = L"J_x = \text{sin}θ", fontsize = 24)
#text!(axxxloc,.7,0.55; text = L"J_z = \text{cos}θ", fontsize = 24)

axislegend(axxxloc)
#axislegend(axxxnonloc, position = :rb)
#axislegend(axzzloc)
#axislegend(axzznonloc, position = :rb)
display("image/png",fig)


## determine T at which data shows 1% deviation on average

Tdifarr = []
for j in eachindex(nthetavals)
    t = pi * nthetavals[j]

    PMFRGvals = chixx0_theta[1,:,j]
    dif = [abs.(chixxloc(t,Tvals[i]) - PMFRGvals[i])/PMFRGvals[i] for i in eachindex(Tvals)]
    posdif = findfirst(x->x > 0.05, dif)
    #println("pos = ",posdif)
    #println("dif[posdif] = ", dif[posdif])
    #println("Tvals[posdif] = ", Tvals[posdif])
    push!(Tdifarr, Tvals[posdif])

    if t != 0.0
        PMFRGvals = chixx0_theta[2,:,j]
        dif = [abs.(chixxnonloc(t,Tvals[i]) - PMFRGvals[i])/abs(PMFRGvals[i]) for i in eachindex(Tvals)]
        posdif = findfirst(x->x > 0.05, dif)
        #println("pos = ",posdif)
        #println("dif[posdif] = ", dif[posdif])
        #println("Tvals[posdif] = ", Tvals[posdif])
        push!(Tdifarr, Tvals[posdif])
    end
    
    PMFRGvals = chizz0_theta[1,:,j]
    dif = [abs.(chizzloc(t,Tvals[i]) - PMFRGvals[i])/abs(PMFRGvals[i]) for i in eachindex(Tvals)]
    posdif = findfirst(x->x > 0.05, dif)
    #println("pos = ",posdif)
    #println("dif[posdif] = ", dif[posdif])
    #println("Tvals[posdif] = ", Tvals[posdif])
    push!(Tdifarr, Tvals[posdif])

    if t != 0.5pi
        PMFRGvals = chizz0_theta[2,:,j]
        dif = [abs.(chizznonloc(t,Tvals[i]) - PMFRGvals[i])/abs(PMFRGvals[i]) for i in eachindex(Tvals)]
        posdif = findfirst(x->x > 0.05, dif[1:end])
        #println("pos = ",posdif)
        #println("dif[posdif] = ", dif[posdif])
        #println("Tvals[posdif] = ", Tvals[posdif])
        push!(Tdifarr, Tvals[posdif])
    end

end

##

sum(Tdifarr)/(length(Tdifarr)+1)


## equal time onsite correlations

fig = Figure()
ax = Axis(fig[1,1])

for ntheta in [0.25]
    file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/PMFRG_XXZ_Tflow_2Polymer_J1=1.0_N=12_NLen=2_theta=$(strno(ntheta,2))pi.h5")
    Tvals = read(file,"State/Tvals")
    chixx_R_nu_T = read(file,"State/chixx")
    chizz_R_nu_T = read(file,"State/chizz")

    chixx_R_T = [chixx_R_nu_T[i,1,j] + 0*2*sum(chixx_R_nu_T[i,2:end,j]) for i in 1:2, j in eachindex(Tvals)]
    chizz_R_T = [chizz_R_nu_T[i,1,j] + 0*2*sum(chizz_R_nu_T[i,2:end,j]) for i in 1:2, j in eachindex(Tvals)]

    scatter!(ax,Tvals,chixx_R_T[1,:], label = "\text{xx}")
    scatter!(ax,Tvals,chizz_R_T[1,:], label = "\text{zz}")
end

xlims!(ax,0.0,1.0)
display("image/png",fig)


## vertex symmetry u <-> t
using HDF5, CairoMakie

folder = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/"
filenames = readdir(folder)
filenames = filter!(x->x[1]=='v',filenames)

figD = Figure()
axD = Axis(figD[1,1], ylabel = L"V_\text{zxzx}", xlabel = L"T/J")
maxV = 0.0

for j in eachindex(filenames)
    file = h5open(folder*filenames[j])
    Tvals = read(file,"State/Tvals")
    N = read(file,"Par/NumericalParams/N")
    Vzxzx = read(file,"State/Vzxzx")
    Vzzxx = read(file,"State/Vzzxx")

    for is in 1:N
        for it in 1:N
            for iu in 1:N
                for Rij in 1:1
                    Dzxzx = Vzxzx[Rij,is,it,iu,:] + Vzzxx[Rij,it,is,iu,:]
                    if is == N && it == N && iu == N && Rij == 1 && j == 4
                        scatter!(axD,Tvals,Vzxzx[Rij,is,it,iu,:], color = :grey, label = L"V^{T,ii}_\text{zxzx}(s,t,u)")
                        scatter!(axD,Tvals,Dzxzx, color = :black, label = L"V^{T,ii}_\text{zxzx}(s,t,u) + V^{T,ii}_\text{zxzx}(t,s,u)")
                    else
                        scatter!(axD,Tvals,Vzxzx[Rij,is,it,iu,:], color = :grey)
                        scatter!(axD,Tvals,Dzxzx, color = :black)
                    end
                end
            end
        end
    end
end

xlims!(0.0,1.0)
ylims!(axD,-4.0,4.0)
axislegend(axD)

#save("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/plots/vertices/SimpleCubic_Heisenberg_Vzxzx.png",figD)
figD
