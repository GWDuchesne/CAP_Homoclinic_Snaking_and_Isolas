module SnakingAndIsolasProof
# Module for the proofs
# ====================================================================================
# Packages
# ====================================================================================
using RadiiPolynomial, LinearAlgebra, JLD2, Plots, TickTock , Polynomials, ProgressBars, LaTeXStrings
using JuMP
using JuMP: value
using Ipopt
# ====================================================================================
# Module
# ====================================================================================
    include("Other_functions.jl")
    include("Map_functions_Orientable.jl")
    include("Map_derivative_Orientable.jl")
    include("Newton_functions.jl")
    include("Finite_difference.jl")
    include("Map_derivative_Orientable_vector.jl")
    include("Bounds_Functions.jl")
    include("Bound_Final_v2.jl")
# ====================================================================================
# Data
# ====================================================================================
    setdisplay(:midpoint; decorations = false, ng_flag = false)

    export Snaking_Proof
    
    function Snaking_Proof( Eq_name::String, branch_prefix::String , branch_number::Int, ν₁::Float64, ν₂::Float64, ν₃::Float64   )

    dir_name = joinpath(@__DIR__, "Data", Eq_name, "GRID", branch_prefix)
    #dir_name = joinpath(@__DIR__, "Data", "SH Orientable Short", "GRID","Data_SH_Orientable_Grid_")

    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = file["X_grid"]
    η_grid = file["η_grid"]
    N = file["N"]
    κ = file["κ"]
    ξ₁ᵏ = file["ξ₁ᵏ"]
    ξ₂ᵏ = file["ξ₂ᵏ"]
    star = file["star"]
    lengthᵥ = file["lengthᵥ"]
    M₁ = file["M₁"]
    M₂ = file["M₂"]
    νᵧ = interval(ν₁ ) 
    νᵥ = interval(ν₁) 
    ν_w = interval(ν₁)
    νₐ = interval(ν₂) 
    νₚ = interval(ν₃)  

    tick()
    X_cheb , η_cheb, X_grid , η_grid, N, κ , ξ₁ᵏ, ξ₂ᵏ, star, lengthᵥ, M₁, M₂ , Y₀ , Z₀, Z₁ , Z₂r , root_all , νᵧ, νᵥ , νₐ , ν_w ,  νₚ , weight_x =  Bounds_Final_v2!(X_grid,η_grid,N ,κ,ξ₁ᵏ,ξ₂ᵏ,star,lengthᵥ,M₁,M₂, νᵧ, νᵥ, ν_w, νₐ, νₚ )
    tock()

    dir_name = joinpath(@__DIR__, "Data", Eq_name , "PROVEN", branch_prefix)
    mkpath(dir_name)
    jldsave(dir_name*"_proof_$(branch_number).jld2"; X_cheb , η_cheb, X_grid , η_grid, N, κ , ξ₁ᵏ, ξ₂ᵏ, star, lengthᵥ, M₁, M₂ , Y₀ , Z₀, Z₁ , Z₂r , root_all , νᵧ, νᵥ , νₐ , ν_w ,  νₚ , weight_x)
    end
end