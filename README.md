# CAP_Homoclinic_Snaking_and_Isolas
Codes to accompany the work “From heteroclinic loops to homoclinic snaking in reversible systems: rigorous forcing through computer-assisted proofs”

- Solution_Finder_Helper_Gray-Scott.jl: Assists in finding a numerical candidate solution to the Gray-Scott equation.
- PAC.jl: Implements the pseudo-arclength continuation method to generate a data grid suitable for constructing Chebyshev polynomial approximations. 
- Newton_Kantorovich.jl: Rigorously computes the bound Y_0, Z_0, Z_1, Z_2 and the radius r_0 for each segment along a continuation loop.
- r_verification.jl: Verifies that the radii r_0 obtained from Newton_Kantorovich.jl satisfy the radii polynomial condition and are thus valid.
- smoothness_verification.jl: Checks that the smoothness condition holds across all segments.
- phase_verification.jl: Compares the endpoints of the lifted phase.


