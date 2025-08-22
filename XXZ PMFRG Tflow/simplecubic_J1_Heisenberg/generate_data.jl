
using SpinFRGLattices, PMFRG_XXZ_Tflow
using SpinFRGLattices.SimpleCubic

NLenvals = [6]

for NLen in NLenvals
    System = getCubic(NLen,[-sqrt(2),0.0])

    Par = ParamsXXZ(
            System,
            theta = 0.25 * pi,
            T_max = 1e+4,
            T_min = 1e-2,
            N = 8,
            accuracy = 1e-4
        )

    file = "/Users/yannikschaden/Julia/Projects/XXZ PMFRG Tflow/simplecubic_J1_Heisenberg/data/local"*fng(Par)
    SolveFRGXXZ(Par, mainfile = file, dtmin=0.001, Katanin = true);
end
