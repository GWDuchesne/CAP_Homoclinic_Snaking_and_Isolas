# =============================================================
# Verify the condition on the endpoints of the lifted phase
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

function car2polar( y , x)
    if  strictprecedes(-x,interval(0))
        atan.(  y ./ x  )
    elseif strictprecedes(x,interval(0)) && precedes(-y,interval(0)) 
        atan.(  y ./ x  ) .+ interval(π)
    elseif strictprecedes(x,interval(0)) && strictprecedes(y,interval(0)) 
        atan.( y / x ) .- interval(π)
    elseif x == 0 && precedes(-y,interval(0)) 
        interval(π) ./ interval(2) 
    elseif x == 0 && strictprecedes(y,interval(0)) 
        - interval(π) ./ interval(2) 
    else 
        interval(0)
    end
end

logocolors = Colors.JULIA_LOGO_COLORS
colors = [logocolors.blue, logocolors.red, logocolors.green, logocolors.purple]


# SH Snaking 
dir_name = joinpath(@__DIR__, "Data", "SH_Snake", "PROVEN","Data_SH_Snake_Grid__proof_")
printstyled("Verifying θ(0) - θ(1) for SH Snaking\n"; color = :blue)

θ_at_0 = interval(0.)
θ_at_last = interval(0.)
P = plot()
for i = 1:6
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = map( X -> interval.(X),  file["X_grid"])
    X_cheb =   dft_grid2cheb(X_grid)  
    t = range(-1,1, 100 )
    θ₁ = real.( component(X_cheb,16)[1])
    θ₂ = real.( component(X_cheb,17)[1])
    vect_angle = interval.( zeros(100) )
    vect_angle[1] = car2polar( θ₁(interval(-1)) , θ₂(interval(-1)))
    if i == 1
        global θ_at_0 = car2polar( θ₁(interval(-1)) , θ₂(interval(-1)))
    end
    vect_angle[1] = car2polar( θ₁(interval(-1)) , θ₂(interval(-1)))
    for j = 2:100
        vect_angle[j] = car2polar( θ₁(interval(t[j])) , θ₂(interval(t[j])))
        if sup( norm( vect_angle[j] - vect_angle[j-1] ) ) > π
            vect_angle[j] = vect_angle[j] + ExactReal(2)*interval(π)
        end
    end
    plot!(t .+ 2*(i-1) ,mid.(vect_angle), legend = false)
    if i == 6
        global θ_at_last = vect_angle[end]
    end
end

display(P)


sdfsdf
# GS Snaking 

dir_name = joinpath(@__DIR__, "Data", "GS_Snake", "PROVEN","Data_GS_Snake_Grid__proof_")
printstyled("Verifying  r > 0 for GS Snaking\n"; color = :blue)
for i = 1:19
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    X_grid = file["X_grid"]
    X_cheb =   dft_grid2cheb(X_grid)  
    t = range(-1,1, 100 )
    α = mid.( real.( component(X_cheb,1)[1].(t)))
    θ₁ = ( real.( component(X_cheb,16)[1].(t)))
    θ₂ = ( real.( component(X_cheb,17)[1].(t)) )
    r = file["root_all"]
    r₀ = r[1]
    if r₀ > 0 
        printstyled("Branch $(branch_number): Succeed\n"; color = :green)
    else
        printstyled("Branch $(branch_number): Failed\n"; color = :red)
    end
end


# SH Isolas

dir_name = joinpath(@__DIR__, "Data", "SH_Isolas", "PROVEN","Data_SH_Isolas_Grid__proof_")
printstyled("Verifying  r > 0 for SH Isolas\n"; color = :blue)
for i = 1:3
    branch_number = i
    file = jldopen(dir_name*"$(branch_number).jld2")
    r = file["root_all"]
    r₀ = r[1]
    if r₀ > 0 
        printstyled("Branch $(branch_number): Succeed\n"; color = :green)
    else
        printstyled("Branch $(branch_number): Failed\n"; color = :red)
    end
end
