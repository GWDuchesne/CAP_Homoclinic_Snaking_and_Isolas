# =====================================================================================
# Rigorous computation of the bounds Y₀, Z₀, Z₁ and Z₂ and radius r₀ for each segment.
# =====================================================================================

include("SnakingAndIsolasProof.jl")
using .SnakingAndIsolasProof

# SH Snaking
for i =1:6
    Snaking_Proof("SH_Snake","Data_SH_Snake_Grid_",i, 1.4, 1.1, 1.4)
    GC.gc()
end

# GS Snaking 
for i = 1:19
    Snaking_Proof("GS_Snake","Data_GS_Snake_Grid_",i, 1.4, 1.1, 1.4)
    GC.gc()
end

# SH Isolas
for i = 1:19
    if i in [10 ; 15 ;17 ; 18 ; 19] 
        Snaking_Proof("SH_Isolas","Data_SH_Isolas_Grid_",i, 1.5, 1.1, 1.5)
    else
        Snaking_Proof("SH_Isolas","Data_SH_Isolas_Grid_",i, 1.4, 1.1, 1.4)
    end
    GC.gc()
end
