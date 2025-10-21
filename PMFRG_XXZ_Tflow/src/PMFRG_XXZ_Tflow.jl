module PMFRG_XXZ_Tflow
    using RecursiveArrayTools, SpinFRGLattices, HDF5, OrdinaryDiffEq, DiffEqCallbacks, StructArrays
    using SpinFRGLattices.StaticArrays

    export SolveFRGXXZ,ParamsXXZ,OneLoopParamsXXZ,BS3,Vern7,DP5,version,getChixx,getChizz,fng

    include("1_structs.jl")
    include("2_vertices.jl")
    include("3_floweq.jl")
    include("4_susceptibility.jl")
    include("5_savedata.jl")
    include("6_solver.jl")
end # module
