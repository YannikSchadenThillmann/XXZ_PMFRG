
using SpinFRGLattices, PMFRG_XXZ_Tflow
include("/Users/yannikschaden/Julia/Helpers/misc.jl")

System = getPolymer(2,coup = [0.0, 1.0])
nthetavals = 0.25:0.25:0.25

logrange_a(min,max,n,a) = (exp(x)+a for x in LinRange(log(min-a),log(max-a),n))
fine_range = collect(LinRange(0.0,0.05,100))
dense_range = collect(LinRange(0.05,2.,50))
medium_range = collect(logrange_a(2.,10.,20,1.7))
sparse_range = collect(logrange_a(10.,1e+4,10,8.8))
ObsSaveat = unique!(append!(fine_range,dense_range,medium_range,sparse_range))

for ntheta in [-0.027]
    Par = ParamsXXZ(
        System,
        theta = ntheta * pi,
        T_max = 1e+4,
        T_min = 1e-3,
        N = 8,
        accuracy = 1e-3
    )

    file = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/test_"*fng(Par)
    SolveFRGXXZ(Par, mainfile = file, dtmin=0.01, vertexoutput = false);
end

## eval test

using HDF5, CairoMakie

file = h5open("/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/dimer/data/test_ObsSaveat_PMFRG_XXZ_Tflow_2Polymer_J1=1.0_N=8_theta=-0.03pi.h5")
Tvals = read(file,"State/Tvals")
Tvals = filter(x->x<0.05,Tvals)

fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax,Tvals,zeros(length(Tvals)))
fig
