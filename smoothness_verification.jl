# =============================================================
# Verify the smoothness condition
# =============================================================
using RadiiPolynomial, LinearAlgebra, JLD2, Plots , Polynomials
# ====================================================================================
# Module
# ====================================================================================
    include("Other_functions.jl")
    include("Map_functions_Orientable.jl")
    include("Map_derivative_Orientable.jl")
    include("Newton_functions.jl")
    include("Finite_difference.jl")
    include("Map_derivative_Orientable_vector.jl")

# ====================================================================================
# Data
# ====================================================================================
grid2cheb(x_fft::Vector{<:Interval}, N) =
[ ifft!(complex.(getindex.(x_fft, i)), Chebyshev(N)) for i ∈ eachindex(x_fft[1])]

grid2cheb(x_fft::Vector{<:Vector}, N) =
[ifft!(getindex.(x_fft, i), Chebyshev(N)) for i ∈ eachindex(x_fft[1])]

grid2cheb(x_fft::Vector{<:Matrix}, N) =
[ifft!(complex.(getindex.(x_fft, i, j)), Chebyshev(N)) for i ∈ axes(x_fft[1], 1), j ∈ axes(x_fft[1], 2)]

grid2cheb(x_fft::Vector{<:LinearOperator}, N) =
[ifft!(complex.(getindex.(x_fft, i, j)), Chebyshev(N)) for i ∈ indices(codomain(x_fft[1])), j ∈ indices(domain(x_fft[1]))]

grid2cheb(x_fft::Vector{<:Sequence}, N) =
[ifft!(complex.(getindex.(x_fft, i)), Chebyshev(N)) for i ∈ indices(space(x_fft[1]))]   

function cheb2grid(x, N_fft)
vals = fft.(x, N_fft)
return [getindex.(vals, i) for i ∈ eachindex(first(vals))]
end
# =====================================================================

        
function idft!(v_grid::Vector)
    N = length(v_grid) - 1
    N_dft = 2*N
    c = zeros(eltype(v_grid), Chebyshev(N))
    for n ∈ 0:N
        c[n] = v_grid[N+1] + (-1)^n * v_grid[1]
        for k ∈ 1:N-1
            c[n] += 2 * v_grid[N+1-k] * cospi(k*n/N)
        end
        c[n] /= N_dft
    end
    c[N] /= 2
    return c
end
function idft!(v_grid::Vector{Complex{Interval{T}}}) where {T}
    N = length(v_grid) - 1
    N_dft = 2*N
    c = zeros(Complex{Interval{T}}, Chebyshev(N))
    for n ∈ 0:N
        c[n] = v_grid[N+1] + ExactReal((-1)^n) * v_grid[1]
        for k ∈ 1:N-1
            c[n] += ExactReal(2) * v_grid[N+1-k] * setprecision(BigFloat, 256) do
                cospi(interval(BigFloat, k*n // N))
            end
        end
        c[n] /= ExactReal(N_dft)
    end
    c[N] /= ExactReal(2)
    return c
end
dft_grid2cheb(A_grid::Vector{<:LinearOperator}) =
LinearOperator(domain(A_grid[1]), codomain(A_grid[1]),
    [idft!(getindex.(A_grid, i, j)) for i ∈ indices(codomain(A_grid[1])), j ∈ indices(domain(A_grid[1]))])


dft_grid2cheb(x_grid::Vector{<:Sequence}) =
Sequence(space(x_grid[1]), [idft!(getindex.(x_grid, i)) for i ∈ indices(space(x_grid[1]))])


dft_grid2cheb(A_grid::Vector{<:Matrix}) =
[idft!(getindex.(A_grid, i, j)) for i ∈ axes(A_grid[1], 1), j ∈ axes(A_grid[1], 2)]

logocolors = Colors.JULIA_LOGO_COLORS
colors = [logocolors.blue, logocolors.red, logocolors.green, logocolors.purple]
# SH Snaking 
dir_name = joinpath(@__DIR__, "Data", "SH_Snake", "PROVEN","Data_SH_Snake_Grid__proof_")
printstyled("Verifying Smothness for SH Snaking\n"; color = :blue)
for i = 1:6
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = file["X_grid"]
    η_grid = file["η_grid"]
    r = file["root_all"]
    r₀ = interval(r[1])
    X_cheb = dft_grid2cheb(map(X -> interval.(X), X_grid))
    η_cheb = dft_grid2cheb(map(X -> interval.(X), η_grid))
    DX_cheb = differentiate.(X_cheb)
    Dη_cheb = differentiate.(η_cheb)
    f_s = sum( DX_cheb .* η_cheb )
    s =  range(-1,1,1000)
    norm_f = real(f_s(interval( s[1] , s[2])) )
    for j = 3:1000
        s_int = interval( s[j-1] , s[j])
        norm_f = min(norm_f, real(f_s(s_int)) )
    end
    if sup( - norm_f +  r₀ * sum( norm.( Dη_cheb , 1 ))) < 0 
        printstyled("Branch $(branch_number): Succeed\n"; color = :green)
    else
        printstyled("Branch $(branch_number): Failed\n"; color = :red)
    end

end

# GS Snaking 
dir_name = joinpath(@__DIR__, "Data", "GS_Snake", "PROVEN","Data_GS_Snake_Grid__proof_")
printstyled("Verifying Smothness for GS Snaking\n"; color = :blue)
for i = 1:19
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = file["X_grid"]
    η_grid = file["η_grid"]
    r = file["root_all"]
    r₀ = interval(r[1])
    X_cheb = dft_grid2cheb(map(X -> interval.(X), X_grid))
    η_cheb = dft_grid2cheb(map(X -> interval.(X), η_grid))
    DX_cheb = differentiate.(X_cheb)
    Dη_cheb = differentiate.(η_cheb)
    f_s = sum( DX_cheb .* η_cheb )
    s =  range(-1,1,1000)
    norm_f = real(f_s(interval( s[1] , s[2])) )
    for j = 3:1000
        s_int = interval( s[j-1] , s[j])
        norm_f = min(norm_f, real(f_s(s_int)) )
    end
    if sup( - norm_f +  r₀ * sum( norm.( Dη_cheb , 1 ))) < 0 
        printstyled("Branch $(branch_number): Succeed\n"; color = :green)
    else
        printstyled("Branch $(branch_number): Failed\n"; color = :red)
    end
end

# SH Isolas

dir_name = joinpath(@__DIR__, "Data", "SH_Isolas", "PROVEN","Data_SH_Isolas_Grid__proof_")
printstyled("Verifying Smothness for SH Isolas\n"; color = :blue)
for i = 1:19
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = file["X_grid"]
    η_grid = file["η_grid"]
    r = file["root_all"]
    r₀ = interval(r[1])
    X_cheb = dft_grid2cheb(map(X -> interval.(X), X_grid))
    η_cheb = dft_grid2cheb(map(X -> interval.(X), η_grid))
    DX_cheb = differentiate.(X_cheb)
    Dη_cheb = differentiate.(η_cheb)
    f_s = sum( DX_cheb .* η_cheb )
    s =  range(-1,1,1000)
    norm_f = real(f_s(interval( s[1] , s[2])) )
    for j = 3:1000
        s_int = interval( s[j-1] , s[j])
        norm_f = min(norm_f, real(f_s(s_int)) )
    end
    if sup( - norm_f +  r₀ * sum( norm.( Dη_cheb , 1 ))) < 0 
        printstyled("Branch $(branch_number): Succeed\n"; color = :green)
    else
        printstyled("Branch $(branch_number): Failed\n"; color = :red)
    end
end
