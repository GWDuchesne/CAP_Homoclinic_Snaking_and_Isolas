# ===============================================================
# Function computing the bounds Y₀, Z₀, Z₁ and Z₂ and radius r₀
# ===============================================================
function Bounds_Final_v2!(X_grid::Vector{<:Sequence},η_grid::Vector{<:Sequence},N::Int64 ,κ::Int64,ξ₁ᵏ::ComplexF64,ξ₂ᵏ::ComplexF64,star::Union{Vector,Int64},lengthᵥ::Float64,M₁::Array,M₂::Array, νᵧ::Interval{Float64}, νᵥ::Interval{Float64}, ν_w::Interval{Float64}, νₐ::Interval{Float64}, νₚ::Interval{Float64}  )

# ====================================================================================
# norms
# ====================================================================================
α = real(component(X_grid[1],1))[1]
γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(component(X_grid[1],2:29))
η = η_grid[1]

𝑿 = NormedCartesianSpace( ( ℓ∞(), 
                            ℓ¹(GeometricWeight(νᵧ)), ℓ¹(GeometricWeight(νᵧ)),ℓ¹(GeometricWeight(νᵧ)),ℓ¹(GeometricWeight(νᵧ)),
                            ℓ∞(), ℓ¹(GeometricWeight(νᵥ)), ℓ¹(GeometricWeight(νᵥ)),ℓ¹(GeometricWeight(νᵥ)),ℓ¹(GeometricWeight(νᵥ)),
                            ℓ¹(GeometricWeight(ν_w)), ℓ¹(GeometricWeight(ν_w)),ℓ¹(GeometricWeight(ν_w)),ℓ¹(GeometricWeight(ν_w)),
                            ℓ∞(), ℓ∞(),ℓ∞(), ℓ∞(), ℓ∞(),
                            ℓ¹(GeometricWeight(νₐ)), ℓ¹(GeometricWeight(νₐ)),ℓ¹(GeometricWeight(νₐ)),ℓ¹(GeometricWeight(νₐ)),
                            ℓ∞(), 
                            ℓ∞(), 
                            ℓ¹(GeometricWeight(νₚ)), ℓ¹(GeometricWeight(νₚ)),ℓ¹(GeometricWeight(νₚ)),ℓ¹(GeometricWeight(νₚ)),
                            )  , ℓ∞())

𝒀 = NormedCartesianSpace( ( ℓ∞(), 
                            ℓ¹(GeometricWeight(νᵧ)), ℓ¹(GeometricWeight(νᵧ)),ℓ¹(GeometricWeight(νᵧ)),ℓ¹(GeometricWeight(νᵧ)),
                            ℓ∞(), ℓ¹(GeometricWeight(νᵥ)), ℓ¹(GeometricWeight(νᵥ)),ℓ¹(GeometricWeight(νᵥ)),ℓ¹(GeometricWeight(νᵥ)),
                            ℓ¹(GeometricWeight(ν_w)), ℓ¹(GeometricWeight(ν_w)),ℓ¹(GeometricWeight(ν_w)),ℓ¹(GeometricWeight(ν_w)),
                            ℓ∞(), ℓ∞(),ℓ∞(), ℓ∞(),
                            ℓ¹(GeometricWeight(νₐ)), ℓ¹(GeometricWeight(νₐ)),ℓ¹(GeometricWeight(νₐ)),ℓ¹(GeometricWeight(νₐ)),
                            ℓ∞(),  ℓ∞(), 
                            ℓ∞(), 
                            ℓ¹(GeometricWeight(νₚ)), ℓ¹(GeometricWeight(νₚ)),ℓ¹(GeometricWeight(νₚ)),ℓ¹(GeometricWeight(νₚ)),
                            )  , ℓ∞())

𝑿ᵧ = ℓ¹(GeometricWeight(νᵧ))
𝑿ᵥ = ℓ¹(GeometricWeight(νᵥ))
𝑿_w = ℓ¹(GeometricWeight(ν_w))
𝑿ₚ = ℓ¹(GeometricWeight(νₚ))
𝑿ₐ = ℓ¹(GeometricWeight(νₐ))

𝑿ᵧ⁴ = NormedCartesianSpace( 𝑿ᵧ , ℓ∞() ) 

space_γ_int = CosFourier( order(γ)[1] , interval(1) ) × SinFourier( order(γ)[2] , interval(1) ) × CosFourier( order(γ)[3] , interval(1) ) × SinFourier( order(γ)[4] , interval(1) ) 
space_v_int = Fourier( order(v)[1] , interval(1) ) × Fourier( order(v)[2] , interval(1) ) × Fourier( order(v)[3] , interval(1) ) × Fourier( order(v)[4] , interval(1) ) 
space_w_int = (Fourier( order(w)[1][1] , interval(1) ) ⊗ Taylor( order(w)[1][2]  ) ) ×   (Fourier( order(w)[2][1] , interval(1) ) ⊗ Taylor( order(w)[2][2]  ) ) ×   (Fourier( order(w)[3][1] , interval(1) ) ⊗ Taylor( order(w)[3][2]  ) ) ×   (Fourier( order(w)[4][1] , interval(1) ) ⊗ Taylor( order(w)[4][2]  ) ) 
space_x_int = space_γ_int × ParameterSpace() × space_v_int × space_w_int × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace() ×  ParameterSpace()^10 × space(p)

α_int = interval(α)
γ_int = Sequence(space_γ_int,interval.(γ[:]))
λ_int = interval(λ)
v_int = Sequence(space_v_int,interval.(v[:]))
w_int = Sequence(space_w_int,interval.(w[:]))
L_int = interval(L)
θ₀_int = interval.(θ₀)
σ_int = interval.(σ)
a_int = interval.(a)
ω_int = interval(ω)
eigenpairs_int = interval.(eigenpairs)
p_int = interval.(p)
x_int = Sequence(space_x_int, [ γ_int[:]; λ_int; v_int[:]; w_int[:]; L_int; θ₀_int[:]; σ_int[:] ; a_int[:]; ω_int; eigenpairs_int[:]; p_int[:]])

space_f_int = space_f!(x_int)
space_X_int = ParameterSpace() × space_x_int
X_int = Sequence( space_X_int, [ α_int ; x_int[:] ] )
η_int = Sequence( space_X_int,  interval.(η[:])  )

space_F_int = space_F!(X_int)

d₁ = [0]
for n  in axes(M₁,1)
    for m  in axes(M₁,2)
        for ℓ  in axes(M₁,3)
            if M₁[n,m,ℓ] != 0
                d₁[1] = maximum([d₁;(n-1)+m-1])
            end
        end
    end
end
d₁ = d₁[1]
d₂ = [0]
for n  in axes(M₂,1)
    for m  in axes(M₂,2)
        for ℓ  in axes(M₂,3)
            if M₂[n,m,ℓ] != 0
                d₂[1] = maximum([d₂;(n-1)+m-1])
            end
        end
    end
end
d₂ = d₂[1]
space_Γ_int_ext = SinFourier( order(γ)[1] , interval(1) ) × CosFourier( d₁*order(γ)[2] , interval(1) ) × SinFourier( order(γ)[3] , interval(1) ) × CosFourier( d₂*order(γ)[4] , interval(1) ) 
space_V_int_ext = ParameterSpace() × Fourier( order(v)[1] , interval(1) ) × Fourier( d₁*order(v)[2] , interval(1) ) × Fourier( order(v)[3] , interval(1) ) × Fourier( d₂*order(v)[4] , interval(1) ) 
space_W_int_ext = (Fourier( order(w)[1][1] , interval(1) ) ⊗ Taylor( order(w)[1][2]  ) ) ×   (Fourier( d₁*order(w)[2][1] , interval(1) ) ⊗ Taylor( d₁*order(w)[2][2]  ) ) ×   (Fourier( order(w)[3][1] , interval(1) ) ⊗ Taylor( order(w)[3][2]  ) ) ×   (Fourier( d₂*order(w)[4][1] , interval(1) ) ⊗ Taylor( d₂*order(w)[4][2]  ) ) 
space_G_int_ext = ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × Chebyshev( order(a)[1] + 1 ) × Chebyshev( d₁*order(a)[2] + 1 )  × Chebyshev( order(a)[3] + 1 ) × Chebyshev( d₂*order(a)[4] + 1 )   × ParameterSpace() × ParameterSpace()
space_E_int_ext = ParameterSpace()^10
space_P_int_ext = (Taylor( order(p)[1][1]) ⊗ Taylor( order(p)[1][2])) × (Taylor( d₁*order(p)[2][1]) ⊗ Taylor( d₁*order(p)[2][2])) × (Taylor( order(p)[3][1]) ⊗ Taylor( order(p)[3][2])) × (Taylor( d₂*order(p)[4][1]) ⊗ Taylor( d₂*order(p)[4][2])) 
space_F_int_ext = ParameterSpace() × space_Γ_int_ext × space_V_int_ext × space_W_int_ext × space_G_int_ext × space_E_int_ext × space_P_int_ext

M₁_int = interval.(M₁)
M₂_int = interval.(M₂)
ξ₁ᵏ_int = interval(ξ₁ᵏ)
ξ₂ᵏ_int = interval(ξ₂ᵏ)
lengthᵥ_int = interval(lengthᵥ)

d = max(d₁,d₂)
space_γ_int_ext = CosFourier( d*order(γ)[1] , interval(1) ) × SinFourier( order(γ)[2] , interval(1) ) × CosFourier( d*order(γ)[3] , interval(1) ) × SinFourier(order(γ)[4] , interval(1) ) 
space_v_int_ext = Fourier( d*order(v)[1] , interval(1) ) × Fourier( order(v)[2] , interval(1) ) × Fourier( d*order(v)[3] , interval(1) ) × Fourier( order(v)[4] , interval(1) ) 
space_w_int_ext = (Fourier( d*order(w)[1][1] , interval(1) ) ⊗ Taylor( d*order(w)[1][2]  ) ) ×   (Fourier( order(w)[2][1] , interval(1) ) ⊗ Taylor( order(w)[2][2]  ) ) ×   (Fourier( d*order(w)[3][1] , interval(1) ) ⊗ Taylor( d*order(w)[3][2]  ) ) ×   (Fourier( order(w)[4][1] , interval(1) ) ⊗ Taylor( order(w)[4][2]  ) ) 
space_a_int_ext = Chebyshev(d*order(a)[1]) × Chebyshev(order(a)[2]) × Chebyshev(d*order(a)[3]) × Chebyshev(order(a)[4]) 
space_p_int_ext = (Taylor( d*order(p)[1][1]) ⊗ Taylor( d*order(p)[1][2])) × (Taylor( order(p)[2][1]) ⊗ Taylor( order(p)[2][2])) × (Taylor( d*order(p)[3][1]) ⊗ Taylor( d*order(p)[3][2])) × (Taylor( order(p)[4][1]) ⊗ Taylor( order(p)[4][2])) 

space_X_int_ext = ParameterSpace() × space_γ_int_ext × ParameterSpace() × space_v_int_ext × space_w_int_ext × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space_a_int_ext × ParameterSpace() ×  ParameterSpace()^10 × space_p_int_ext

space_γ_ext = CosFourier( d*order(γ)[1] , 1 ) × SinFourier( order(γ)[2] , 1 ) × CosFourier( d*order(γ)[3] , 1) × SinFourier(order(γ)[4] , 1 ) 
space_v_ext = Fourier( d*order(v)[1] , 1 ) × Fourier( order(v)[2] , 1 ) × Fourier( d*order(v)[3] , 1 ) × Fourier( order(v)[4] , 1 ) 
space_w_ext = (Fourier( d*order(w)[1][1] , 1 ) ⊗ Taylor( d*order(w)[1][2]  ) ) ×   (Fourier( order(w)[2][1] , 1 ) ⊗ Taylor( order(w)[2][2]  ) ) ×   (Fourier( d*order(w)[3][1] , 1 ) ⊗ Taylor( d*order(w)[3][2]  ) ) ×   (Fourier( order(w)[4][1] ,1) ⊗ Taylor( order(w)[4][2]  ) ) 
space_a_ext = Chebyshev(d*order(a)[1]) × Chebyshev(order(a)[2]) × Chebyshev(d*order(a)[3]) × Chebyshev(order(a)[4]) 
space_p_ext = (Taylor( d*order(p)[1][1]) ⊗ Taylor( d*order(p)[1][2])) × (Taylor( order(p)[2][1]) ⊗ Taylor( order(p)[2][2])) × (Taylor( d*order(p)[3][1]) ⊗ Taylor( d*order(p)[3][2])) × (Taylor( order(p)[4][1]) ⊗ Taylor( order(p)[4][2])) 

space_X_ext = ParameterSpace() × space_γ_ext × ParameterSpace() × space_v_ext × space_w_ext × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space_a_ext × ParameterSpace() ×  ParameterSpace()^10 × space_p_ext

# ====================================================================================
# Functions
# ====================================================================================

grid2cheb(x_fft::Vector{<:Interval}, N) =
[rifft!(complex.(getindex.(x_fft, i)), Chebyshev(N)) for i ∈ eachindex(x_fft[1])]

grid2cheb(x_fft::Vector{<:Vector}, N) =
    [ifft!(complex.(getindex.(x_fft, i)), Chebyshev(N)) for i ∈ eachindex(x_fft[1])]

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

# ====================================================================================
# Variables
# ====================================================================================
clc()
printstyled("============================================================\n"; color = :blue)
printstyled("Set up and Converting Grid to Chebyshev\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)
println("▪ Grid to Chebyshev in DFT")
N_fft = nextpow(2, 2N + 1)
npts = N_fft÷2 + 1
space_X = space(X_grid[1])
space_F = space_F!(X_grid[1])

# More functions
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


X_cheb_sequence = dft_grid2cheb(map(X -> interval.(X),X_grid)) 
X_fft_mid = map(X -> mid.(X), cheb2grid(X_cheb_sequence, N_fft))
X_cheb = X_cheb_sequence[:]
X_cheb_sequence = Sequence(space_X_int,X_cheb)

# Checking if ν is not too small for the proof 
#X_cheb_sequence = Sequence( space_X_int, X_cheb)

σ₁_cheb = component(X_cheb_sequence , 18)[1]
σ₂_cheb = component(X_cheb_sequence , 19)[1]
norm_z = sqrt( norm(σ₁_cheb^2 + σ₂_cheb^2,interval(1)) )

θ₁_cheb = component(X_cheb_sequence , 16)[1]
θ₂_cheb = component(X_cheb_sequence , 17)[1]
norm_y = sqrt( norm(θ₁_cheb^2 + θ₂_cheb^2,interval(1)) )

if inf(norm_z) > sup( νₚ ) && inf(norm_y) > sup( ν_w )
    error(String("We need νₚ > $(inf(norm_z)) and ν_w > $(inf(norm_y))" ))
elseif inf(norm_z) > sup( νₚ )
    error(String(" We need νₚ > $(inf(norm_z))"))
elseif inf(norm_y) > sup( ν_w )
    error(String("We need ν_w > $(inf(norm_y))"))
end

η_cheb_sequence = dft_grid2cheb(map(X -> interval.(X),η_grid)) 
η_fft = map(X -> mid.(X), cheb2grid(η_cheb_sequence, N_fft))
η_cheb = η_cheb_sequence[:]

A_fft = A_fft!(  Vector{Vector}( undef, 29 ), space_X, space_F , X_fft_mid, κ, star, lengthᵥ, M₁, M₂ , η_fft ,N_fft)
A_cheb = A_cheb!(Vector{Any}( undef, 29 ) , A_fft , N)
A_fft = nothing
dN = dN!( N , M₁ , M₂ , order(p) , order(w))
dN_ftt = nextpow.(2, 2dN  .+1)

# Some other Variables needed later
Λ_tails = Vector{Vector}(undef, 29)
λ_grid = Vector{Sequence}(undef, lastindex(X_grid))
λ₁_grid = Vector{Sequence}(undef, lastindex(X_grid))
λ₂_grid = Vector{Sequence}(undef, lastindex(X_grid))
for i = 1:lastindex(X_grid)
    λ_grid[i] = component( X_grid[i] , 6)
    λ₁_grid[i] = Sequence( ParameterSpace(),  [component( X_grid[i] , 25)[1]])
    λ₂_grid[i] = Sequence( ParameterSpace(),  [component( X_grid[i] , 25)[6]])
end

λ_cheb = dft_grid2cheb(map(X -> interval.(X),λ_grid)) 
λ₁_cheb = dft_grid2cheb(map(X -> interval.(X),λ₁_grid)) 
λ₂_cheb = dft_grid2cheb(map(X -> interval.(X),λ₂_grid)) 

# ======================================================================================================
#  Y₀  
# ======================================================================================================
# Y₀ body  
printstyled("============================================================\n"; color = :blue)
printstyled("1. Bound Y₀\n"; color = :blue)   
printstyled("============================================================\n"; color = :blue)

println("▪ Computing πᴺAF")
function AF!( A_cheb::Vector, space_X_int::CartesianProduct, space_F_int::CartesianProduct , dN::Vector{Int64},  dN_ftt::Vector{Int64}, X_cheb::Vector{<:Sequence}, η_cheb::Vector{<:Sequence}, κ::Int64, ξ₁ᵏ_int::Complex{Interval{Float64}} ,ξ₂ᵏ_int::Complex{Interval{Float64}}, star::Union{Vector,Int64}, lengthᵥ_int::Interval{Float64}, M₁_int::Array{<:Interval}, M₂_int::Array{<:Interval} )
    F_fft = F_fft!(Vector{Any}( undef, 29 ), space_X_int, space_F_int , dN_ftt, X_cheb, η_cheb, κ, ξ₁ᵏ_int ,ξ₂ᵏ_int, star, lengthᵥ_int, M₁_int, M₂_int)
    AF_fft = AF_fft!(Vector{Vector}( undef, 29 ), A_cheb, F_fft, space_X_int, space_F_int, dN_ftt)
    AF = Sequence( space_X_int , AF_cheb!( AF_fft, dN) )
    return AF
end
AF = AF!( A_cheb, space_X_int, space_F_int , dN,  dN_ftt, X_cheb, η_cheb, κ, ξ₁ᵏ_int ,ξ₂ᵏ_int, star, lengthᵥ_int, M₁_int, M₂_int )
GC.gc()

# Y₀ tails
println(" ")
println("▪ Computing π ᪲ AF")
Y₀_tails_cheb = Y₀_tails_fft!( F_ext_fft!( Vector{Any}( undef, 29 ), X_cheb, η_cheb, space_X_int, space_F_int_ext , dN_ftt, κ, ξ₁ᵏ_int ,ξ₂ᵏ_int, star, lengthᵥ_int, M₁_int, M₂_int) , λ_cheb,λ₁_cheb,λ₂_cheb, dN, dN_ftt,space_F ,w ,a,p,space_F_int_ext, Λ_tails) 
GC.gc()

Y₀_vec = Vector{Any}(undef, 29)
for i = 1:29
    if i in [1;6;15;16;17;18;19;24;25]
        local space_test_in = ℓ∞()
    elseif i in [ 2 ;3;4;5]
        local space_test_in =  ℓ¹(GeometricWeight(νᵧ))
    elseif i in [ 7;8;9;10]
        local space_test_in =  ℓ¹(GeometricWeight(νᵥ))
    elseif i in [ 11;12;13;14]
        local space_test_in =  ℓ¹(GeometricWeight(ν_w))
    elseif i in [ 20;21;22;23]
        local space_test_in =  ℓ¹(GeometricWeight(νₐ))
    elseif i in [ 26;27;28;29]
        local space_test_in =  ℓ¹(GeometricWeight(νₚ))
    end
    if i in [1;6;15;16;17;18;23;24;25]
        local space_test_out = ℓ∞()
    elseif i in [ 2 ;3;4;5]
        local space_test_out =  ℓ¹(GeometricWeight(νᵧ))
    elseif i in [ 7;8;9;10]
        local space_test_out =  ℓ¹(GeometricWeight(νᵥ))
    elseif i in [ 11;12;13;14]
        local space_test_out =  ℓ¹(GeometricWeight(ν_w))
    elseif i in [ 19;20;21;22]
        local space_test_out =  ℓ¹(GeometricWeight(νₐ))
    elseif i in [ 26;27;28;29]
        local space_test_out =  ℓ¹(GeometricWeight(νₚ))
    end
    global Y₀_vec[i] = norm(  component(AF , i ) , space_test_in  ) +  norm( component( Y₀_tails_cheb ,i ) , space_test_out)   
end
AF = nothing
Y₀_tails_cheb = nothing
GC.gc()
# ======================================================================================================
#  Z₁
# ======================================================================================================
printstyled("============================================================\n"; color = :blue)
printstyled("2. Bound Z₁\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)
Aᵢⱼ_cheb = Aᵢⱼ_cheb!(Matrix{Any}( undef, 29 , 29), space_X , space_F, N, N_fft , X_fft_mid, κ, star, lengthᵥ, M₁, M₂ , η_fft)
GC.gc()

# Z₁ tails
println("▪ Computing π ᪲ A[DF(̄x) - A†]")

DN = DN!( N , M₁ , M₂ , order(p) , order(w), κ)
DN_ftt = copy(DN)
for i = 1:29
    for j = 1:29
    if DN[i,j] >0
        DN_ftt[i,j] = nextpow.(2, 2DN[i,j] )
    else
        DN_ftt[i,j] = 1
    end
    end
end
Λ =  Vector{Any}(undef, 29)
for i = 1:29
    Λ[i] = maximum(interval.(Λ_tails[i]))
end
Z₁_₁  =  Z₁_mat!(zeros( Interval{Float64} ,29,29), X_cheb, DN, DN_ftt, Λ , space_X_int, M₁_int, M₂_int , 𝑿ᵧ , 𝑿ᵥ , 𝑿_w , 𝑿ₚ , 𝑿ₐ  ,νₚ , νₐ ,ν_w , Aᵢⱼ_cheb ,space_X_int_ext,space_F_int, νᵥ,νᵧ)
GC.gc()

println( maximum( Z₁_₁ ))


# Z₁ body
println(" ")
println("▪ Computing πᴺA[DF(̄x) - A†]")
Z₁_ = Z₁_body!(copy(Z₁_₁), X_cheb, Aᵢⱼ_cheb, DN, DN_ftt, M₁_int, M₂_int , νᵧ, νᵥ,  νₐ, space_F,space_X_int, ν_w, νₚ)
GC.gc()

# Weights
println(" ")
weight_x = weight_x!(Z₁_ )

# Problem dimensions
m = 29  # number of rows in A
n = 29  # number of variables

# Example matrix A (strictly positive to avoid division by zero)


# Lower bound α for the maximum value constraint
t_min = 0.6  # must satisfy 0 < α < 1

# Define optimization model
model = Model(Ipopt.Optimizer)

# Decision variables:
# v: vector of weights (positive and sum to 1)
# t_max: upper bound on the expressions v[i] * sum_j A[i,j] / v[j]
@variable(model, v[1:n] >= 1e-9)
@variable(model, t_min <= t_max <= 0.999)

# Normalization constraint: sum of v equals 1
@constraint(model, sum(v) == 1)

# For each row i, enforce:
# v[i] * sum_j (A[i,j] / v[j]) ≤ t_max
@NLconstraint(model, [i = 1:m], v[i] * sum(mid.(Z₁_)[i,j] / v[j] for j in 1:n) <= t_max)

# Objective: minimize the maximum value across all rows
@objective(model, Min, t_max)

# Solve the problem
optimize!(model);

# Extract and display results

    v_opt = value.(v)
    t_opt = value(t_max)
    norm(1 ./ v_opt) .* v_opt

   weight_x =  interval.( norm(1 ./ v_opt) .* v_opt)



Y₀ = maximum( Y₀_vec .* weight_x)
printstyled("✶ Y₀ = $(sup(Y₀))\n"; color = :green)

# Z₁ computing with weights  
Z₁ = maximum( (Z₁_ * ( interval(1) ./ weight_x )) .* (   weight_x) )
printstyled("✶ Z₁ = $(sup(Z₁))\n"; color = :green)
GC.gc()
# Z₀ 
printstyled("============================================================\n"; color = :blue)
printstyled("3. Bound Z₀\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)


function space_operatorᵢⱼ!( space_ind, i )
    if space_ind == 1
        if i in [1;6;15;16;17;18;19;24;25]
            space_out = ℓ∞()
        elseif i in [ 2 ;3;4;5]
            space_out =  ℓ¹(GeometricWeight(νᵧ))
        elseif i in [ 7;8;9;10]
            space_out =  ℓ¹(GeometricWeight(νᵥ))
        elseif i in [ 11;12;13;14]
            space_out =  ℓ¹(GeometricWeight(ν_w))
        elseif i in [ 20;21;22;23]
            space_out =  ℓ¹(GeometricWeight(νₐ))
        elseif i in [ 26;27;28;29]
            space_out =  ℓ¹(GeometricWeight(νₚ))
        end
    else
        if i in [1;6;15;16;17;18;23;24;25]
            space_out = ℓ∞()
        elseif i in [ 2 ;3;4;5]
            space_out =  ℓ¹(GeometricWeight(νᵧ))
        elseif i in [ 7;8;9;10]
            space_out =  ℓ¹(GeometricWeight(νᵥ))
        elseif i in [ 11;12;13;14]
            space_out =  ℓ¹(GeometricWeight(ν_w))
        elseif i in [ 19;20;21;22]
            space_out =  ℓ¹(GeometricWeight(νₐ))
        elseif i in [ 26;27;28;29]
            space_out =  ℓ¹(GeometricWeight(νₚ))
        end
    end
    return space_out
end



println("▪ Computing DF on Grid")
DF = DF_ext_fft!(Matrix{Any}(undef,29,29), X_cheb, η_cheb, DN_ftt, space_X_int, space_X_int ,space_F_int, κ, star, lengthᵥ_int, M₁_int, M₂_int )
F = zeros(space_F);

A_ = A_!(zeros( Interval{Float64}, 29 , 29 ) ,Aᵢⱼ_cheb,space_F_int,space_X_int,νᵧ,νᵥ,ν_w,νₐ,νₚ);

DF_cheb_truncated = zeros( Sequence{Chebyshev, Vector{Complex{Interval{Float64}}}} , space_X_int, space_F_int)
DF_tail_norm = zeros( Interval{Float64} , 29,29)
    for i in 1:29
        for j = 1:29
            if j in [1;6;15;16;17;18;19;24;25]
                local space_test_in2 = ℓ∞()
            elseif j in [ 2 ;3;4;5]
                local space_test_in2 =  ℓ¹(GeometricWeight(νᵧ))
            elseif j in [ 7;8;9;10]
                local space_test_in2 =  ℓ¹(GeometricWeight(νᵥ))
            elseif j in [ 11;12;13;14]
                local space_test_in2 =  ℓ¹(GeometricWeight(ν_w))
            elseif j in [ 20;21;22;23]
                local space_test_in2 =  ℓ¹(GeometricWeight(νₐ))
            elseif j in [ 26;27;28;29]
                local space_test_in2 =  ℓ¹(GeometricWeight(νₚ))
            end
            if i in [1;6;15;16;17;18;23;24;25]
                local space_test_out2 = ℓ∞()
            elseif i in [ 2 ;3;4;5]
                local space_test_out2 =  ℓ¹(GeometricWeight(νᵧ))
            elseif i in [ 7;8;9;10]
                local space_test_out2 =  ℓ¹(GeometricWeight(νᵥ))
            elseif i in [ 11;12;13;14]
                local space_test_out2 =  ℓ¹(GeometricWeight(ν_w))
            elseif i in [ 19;20;21;22]
                local space_test_out2 =  ℓ¹(GeometricWeight(νₐ))
            elseif i in [ 26;27;28;29]
                local space_test_out2 =  ℓ¹(GeometricWeight(νₚ))
            end
            if isassigned(DF,i,j)
                component( DF_cheb_truncated ,i,j)[:,:] = grid2cheb( map( X -> X[:,:],DF[i , j])  , DN[i,j] )
                tail_compo = LinearOperator(  space_X_int[j], space_F_int[i],  component( DF_cheb_truncated ,i,j)[:,:] )
                component( DF_cheb_truncated ,i,j)[:,:] = map( X -> project( X , Chebyshev(N)) ,component( DF_cheb_truncated ,i,j)[:,:])
                tail_compo[:,:] .= (x -> (x[0:N] .= 0; x)).(tail_compo[:,:])
                DF_tail_norm[i,j] = opnorm( norm.(tail_compo ,1),space_test_in2,space_test_out2)
            end
        end
    end

DF = nothing
GC.gc()

DF_mid = Matrix{Any}(undef,29,29)
ρ = zeros(Interval{Float64}, 29,29)
DF_mid_norm = zeros(Interval{Float64}, 29,29)

A_mid = Matrix{Any}(undef,29,29)
r = zeros(Interval{Float64}, 29,29)
A_mid_norm = zeros(Interval{Float64}, 29,29)

GC.gc()
println("")
println("▪ Computing Error on the Tail")
iter = ProgressBar(1:29)
set_description(iter, "    - ")
for i in iter
    for j = 1:29
        if DN_ftt[i,j] > 1
            DF_mid[i,j] = cheb2grid( component( DF_cheb_truncated ,i,j)[:,:] ,  nextpow( 2, 4*N + 1 )     ) 
            DF_temp = grid2cheb( DF_mid[i,j] , 2N)
            DF_mid_norm[i,j]   = opnorm(LinearOperator( space_X_int[j], space_F_int[i], norm.( map( X -> interval.(X),map(  X -> mid.(X) ,DF_temp ) ) ,1)) , space_operatorᵢⱼ!( 1, j ), space_operatorᵢⱼ!( 2 , i )  )
            ρ[i,j] = opnorm(LinearOperator( space_X_int[j], space_F_int[i], norm.( map( X -> interval.(X),map(  X -> radius.(X) , DF_temp )).*interval(-1,1) ,1)), space_operatorᵢⱼ!( 1, j ), space_operatorᵢⱼ!( 2 , i )  )
            DF_mid[i,j] = map(X -> mid.(X) , DF_mid[i,j])
        end
        A_mid[i,j] = cheb2grid( Aᵢⱼ_cheb[i,j] ,  nextpow( 2, 4*N + 1 )     ) 
        A_temp = grid2cheb( A_mid[i,j] , 2N)
        A_mid_norm[i,j]  = opnorm(LinearOperator( space_F_int[j], space_X_int[i], norm.(  map( X -> interval.(X),  A_temp ),1)) , space_operatorᵢⱼ!( 2, j ), space_operatorᵢⱼ!( 1 , i )  )
        r[i,j] = opnorm(LinearOperator( space_F_int[j], space_X_int[i], norm.(  map( X -> interval.(X),  map(  X -> radius.(X) , A_temp )).*interval(-1,1) ,1)) , space_operatorᵢⱼ!( 2, j ), space_operatorᵢⱼ!( 1 , i )  )
        A_mid[i,j]  =map(X -> mid.(X) , A_mid[i,j])
    end
end 

GC.gc()

B  = zeros(Interval{Float64}, 29,29)
println("")
println("▪ Computing || I - ADFᴺ||")
iter = ProgressBar(1:29)
set_description(iter, "    - ")
I_ADFᴺᵢⱼ  = Vector{Matrix}(undef, nextpow(2, 4N+1))
for i in iter
    for j= 1:29
        for ℓ = 1:nextpow(2, 4N+1)
            I_ADFᴺᵢⱼ[ℓ] = zeros( Complex{Interval{Float64}} , space_X_int[j], space_X_int[i])[:,:]
            if j == i
                I_ADFᴺᵢⱼ[ℓ] = I_ADFᴺᵢⱼ[ℓ] + UniformScaling(interval(1))
            end
        end
       for k = 1:29
            if DN_ftt[k,j] > 1
            I_ADFᴺᵢⱼ = I_ADFᴺᵢⱼ .- map( X -> interval.(X), A_mid[i,k].*DF_mid[k,j]) 
            end
        end
        B[i,j] = opnorm( LinearOperator(  space_X_int[j], space_X_int[i] , norm.( grid2cheb(  I_ADFᴺᵢⱼ ,2N) , 1  ) ) , space_operatorᵢⱼ!( 1, j ) , space_operatorᵢⱼ!( 1, i ) )
    end
end



GC.gc()

function M_ϵ!(n,ϵ)
    return (ExactReal(1) + ExactReal(2)*ϵ + ϵ^2 )*((ExactReal(1) + ϵ)^(interval(n)) )-ExactReal(1)
end         

ϵ = interval(2.22045)*interval(10)^ExactReal(-16)
M_ϵ = M_ϵ!(maximum( [ length(component(X_grid[1],i)) for i = 1:29] ),ϵ)

Z₀ = maximum( ((B + M_ϵ*A_mid_norm*DF_mid_norm + r*DF_mid_norm +  ρ*A_mid_norm + r*ρ + A_* DF_tail_norm)  * ( interval(1) ./ weight_x )) .* (   weight_x) )


# Z₀ = 0
println("")
 printstyled("✶ Z₀ = $(sup(Z₀))\n"; color = :green)



# ======================================================================================================
#  Z₂
# ======================================================================================================
printstyled("============================================================\n"; color = :blue)
printstyled("4. Bound Z₂\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)

X_cheb_norm = Sequence( space_X_int, norm.( X_cheb, 1))
α_int = real(X_cheb_norm[1])
x = component(X_cheb_norm,2:29)
γ_int,λ,v_int,w_int,L,θ₀_int,σ_int,a_int,ω,eigenpairs_int,p_int = x2var!(x)

γ₁ , γ₂ , γ₃ , γ₄ = eachcomponent(γ_int)
v₁ , v₂ , v₃ , v₄ = eachcomponent(v_int)
w₁ , w₂ , w₃ , w₄ = eachcomponent(w_int)
p₁ , p₂ , p₃ , p₄ = eachcomponent(p_int)
a₁ , a₂ , a₃ , a₄ = eachcomponent(a_int)

c₁  =  Polynomial( [ norm( α_int ) , interval(1) / weight_x[1] ] , :r)
c₂  =  Polynomial( [ norm( γ₁ , 𝑿ᵧ ) , interval(1) / weight_x[2] ] , :r)
c₄  =  Polynomial( [ norm( γ₃ , 𝑿ᵧ ) , interval(1) / weight_x[4] ] , :r)
c₂_ext  =  Polynomial( [ norm( γ₁ , 𝑿ᵧ ) , interval(1) / weight_x[2] ] , :r)
c₄_ext  =  Polynomial( [ norm( γ₃ , 𝑿ᵧ ) , interval(1) / weight_x[4] ] , :r)
c₇  =  Polynomial( [ norm( v₁ , 𝑿ᵥ ) , interval(1) / weight_x[7] ] , :r)
c₉  =  Polynomial( [ norm( v₃ , 𝑿ᵥ ) , interval(1) / weight_x[9] ] , :r)
c₁₁  =  Polynomial( [ norm( w₁ , 𝑿_w ) , interval(1) / weight_x[11] ] , :r)
c₁₃  =  Polynomial( [ norm( w₃ , 𝑿_w ) , interval(1) / weight_x[13] ] , :r)
c₁₅  =  Polynomial( [ norm( L ) , interval(1) / weight_x[15] ] , :r)
c₂₀  =  Polynomial( [ norm( a₁ , 𝑿ₐ ) , interval(1) / weight_x[20] ] , :r)
c₂₂  =  Polynomial( [ norm( a₃ , 𝑿ₐ ) , interval(1) / weight_x[22] ] , :r)
c₂₄ =  Polynomial( [ norm( ω_int ) , interval(1) / weight_x[24] ] , :r)
c₂₅_₂ =  Polynomial( [ norm( eigenpairs_int[2] ) , interval(1) / weight_x[25] ] , :r)
c₂₅_₄ =  Polynomial( [ norm( eigenpairs_int[4] ) , interval(1) / weight_x[25] ] , :r)
c₂₅_₇ =  Polynomial( [ norm( eigenpairs_int[7] ) , interval(1) / weight_x[25] ] , :r)
c₂₅_₉ =  Polynomial( [ norm( eigenpairs_int[9] ) , interval(1) / weight_x[25] ] , :r)
c₂₆  =  Polynomial( [ norm( p₁ , 𝑿ₚ ) , interval(1) / weight_x[26] ] , :r)
c₂₈  =  Polynomial( [ norm( p₃ , 𝑿ₚ ) , interval(1) / weight_x[28] ] , :r)

ΔDF = Matrix{ Any }( undef , 29 , 29 )
ΔDF[:,:] .= Polynomial([interval(0)] , :r)

# i = 2
ΔDF[ 2 , 3 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[24] ] , :r)
ΔDF[ 2 , 24] = Polynomial( [ interval(0) ; interval(1)  / weight_x[3] ] , :r)
# i = 4
ΔDF[ 4 , 5 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[24] ] , :r)
ΔDF[ 4 , 24] = Polynomial( [ interval(0) ; interval(1)  / weight_x[5] ] , :r)
# i = 7 
ΔDF[ 7 , 6 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[7]] , :r)
ΔDF[ 7 , 7 ] = Polynomial( [ interval(0) ; interval(1) / weight_x[6] ] , :r)
ΔDF[ 7 , 8 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[24]] , :r)
ΔDF[ 7 , 24 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[8]] , :r)
# i = 9 
ΔDF[ 9 , 6 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[9]] , :r)
ΔDF[ 9 , 9 ] = Polynomial( [ interval(0) ; interval(1) / weight_x[6] ] , :r)
ΔDF[ 9 , 10 ] = Polynomial( [ interval(0) ; interval(1)  / weight_x[24] ] , :r)
ΔDF[ 9 , 24 ] = Polynomial( [ interval(0) ; interval(1) / weight_x[10] ] , :r)


# i = 11
ΔDF[ 11 , 12 ] = Polynomial( [ interval(0) ; interval(1) / weight_x[24] ] , :r)
ΔDF[ 11 , 24] = Polynomial( [ interval(0) ; interval(1) / weight_x[12] ] , :r)

# i = 13
ΔDF[ 13 , 14 ] = Polynomial( [ interval(0) ; interval(1) / weight_x[24] ] , :r)
ΔDF[ 13 , 24] = Polynomial( [ interval(0) ; interval(1) / weight_x[14] ] , :r)

norm_t = sqrt((interval(1)/(weight_x[18]))^2 + (interval(1)/(weight_x[19]))^2)
norm_h = sqrt((interval(1)/(weight_x[16]))^2 + (interval(1)/(weight_x[17]))^2)

ΔDF[ 15 , 26 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_t/νₚ/(νₚ -norm_z  ) ] , :r)
# i = 16
ΔDF[ 16 , 27 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_t/νₚ/(νₚ -norm_z  ) ] , :r)
# i = 17
ΔDF[ 17 , 28 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_t/νₚ/(νₚ -norm_z  ) ] , :r)
# i = 18
ΔDF[ 18 , 29 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_t/νₚ/(νₚ -norm_z  ) ] , :r)

# i = 15
ΔDF[ 15 , 18 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₁ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[26]  ], :r) 
ΔDF[ 15 , 19 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₁ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[26]  ], :r) 
# i = 16
ΔDF[ 16 , 18 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₂ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[27]  ], :r) 
ΔDF[ 16 , 19 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₂ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[27]  ], :r) 
# i = 17
ΔDF[ 17 , 18 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₃ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[28]  ], :r) 
ΔDF[ 17 , 19 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₃ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[28]  ], :r) 
# i = 18
ΔDF[ 18 , 18 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₄ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[29]  ], :r) 
ΔDF[ 18 , 19 ] = Polynomial( [ interval(0) ;  ExactReal(4)*norm( p₄ , 𝑿ₚ )*norm_t/νₚ/(νₚ - norm_z )/log(ExactReal(2)*νₚ/(norm_z + νₚ ))  +  interval(2)/νₚ/weight_x[29]  ], :r) 

norm_T = interval(1)/νₐ + interval(2)νₐ
D₁₀ = Derivative(1,0)
Dw₁ = abs.(D₁₀*w₁)
Dw₂ = abs.(D₁₀*w₂)
Dw₃ = abs.(D₁₀*w₃)
Dw₄ = abs.(D₁₀*w₄)
Dw = [ Dw₁ ; Dw₂ ;  Dw₃ ; Dw₄ ]
norm_y_poly = Polynomial( [ norm_y ]  , :r)
norm_y_plus_r₁ = Polynomial( [ norm_y ; interval(1)/weight_x[16]]  , :r)
norm_y_plus_r₂ = Polynomial( [ norm_y ; interval(1)/weight_x[17]]  , :r)

# i = 19
ΔDF[ 19 , 15 ] = Polynomial( [ interval(0) ; norm_T/interval(2)/weight_x[21] ] , :r )
ΔDF[ 19 , 21 ] = Polynomial( [ interval(0) ; norm_T/interval(2)/weight_x[15] ] , :r )
ΔDF[ 19 , 11 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_h/ν_w/(ν_w -norm_y  ) ] , :r)
ΔDF_₁₉_₁₆ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₁ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[11]  ], :r) 
ΔDF_₁₉_₁₇ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₁ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[11]  ], :r) 

# i = 20
ΔDF[ 20 , 12 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_h/ν_w/(ν_w -norm_y  ) ] , :r)
ΔDF_₂₀_₁₆ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₂ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[12]  ], :r) 
ΔDF_₂₀_₁₇ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₂ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[12]  ], :r) 

# i = 21
ΔDF[ 21 , 15 ] = Polynomial( [ interval(0) ; norm_T/interval(2)/weight_x[23] ] , :r )
ΔDF[ 21 , 23 ] = Polynomial( [ interval(0) ; norm_T/interval(2)/weight_x[15] ] , :r )
ΔDF[ 21 , 13 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_h/ν_w/(ν_w -norm_y  ) ] , :r)
ΔDF_₂₁_₁₆ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₃ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[13]  ], :r) 
ΔDF_₂₁_₁₇ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₃ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[13]  ], :r) 

# i = 22
ΔDF[ 22 , 14 ] = Polynomial( [ interval(0) ; ExactReal(2)*norm_h/ν_w/(ν_w -norm_y  ) ] , :r)
ΔDF_₂₂_₁₆ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₄ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[14]  ], :r) 
ΔDF_₂₂_₁₇ = Polynomial( [ interval(0) ;  ExactReal(8)*norm( w₄ , 𝑿_w )*norm_h/ν_w/(ν_w - norm_y )/log(ExactReal(10)*ν_w/(norm_y + ν_w ))  +  interval(4)/ν_w/weight_x[14]  ], :r) 

# i = 23
ΔDF[ 23 , 18 ] = Polynomial( [interval(0) ; interval(2)/weight_x[18]  ] ,:r)
ΔDF[ 23 , 19 ] = Polynomial( [interval(0) ; interval(2)/weight_x[19]  ] ,:r)
# i = 24
ΔDF[ 24 , 16 ] = Polynomial( [interval(0) ; interval(2)/weight_x[16]  ] ,:r)
ΔDF[ 24 , 17 ] = Polynomial( [interval(0) ; interval(2)/weight_x[17]  ] ,:r)
# i = 25
ΔDF₂₅ = Matrix{Any}(undef, 10 , 29)
ΔDF₂₅ .= Polynomial( interval(0) , :r)

ΔDF₂₅[2,1] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[2,2] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[3,1] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[3,3] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[4,1] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[4,4] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[5,1] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[5,5] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[7,6] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[7,7] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[8,6] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[8,8] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[9,6] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[9,9] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

ΔDF₂₅[10,6] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)
ΔDF₂₅[10,10] = Polynomial( [interval(0) , interval(1)/weight_x[25] ]  , :r)

for ℓ in axes(M₁,3)
    if size(M₁,1) >= 2 
        if ℓ >= 2
            ΔDF₂₅[3,1] = ΔDF₂₅[3,1]+ abs.(M₁_int)[2,1,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₂  
            ΔDF₂₅[8,1] = ΔDF₂₅[8,1]+ abs.(M₁_int)[2,1,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₇ 
        end
        ΔDF₂₅[3,2] = ΔDF₂₅[3,2]+ abs.(M₁_int)[2,1,ℓ]*c₁^(ℓ-1)
        ΔDF₂₅[8,7] = ΔDF₂₅[8,7]+ abs.(M₁_int)[2,1,ℓ]*c₁^(ℓ-1)
    end
    if size(M₁,2) >= 2 
        if ℓ >= 2
            ΔDF₂₅[3,1] = ΔDF₂₅[3,1]+ abs.(M₁_int)[1,2,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₄  
            ΔDF₂₅[8,1] = ΔDF₂₅[8,1]+ abs.(M₁_int)[1,2,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₇  
        end
        ΔDF₂₅[3,4] = ΔDF₂₅[3,4]+ abs.(M₁_int)[1,2,ℓ]*c₁^(ℓ-1)
        ΔDF₂₅[8,9] = ΔDF₂₅[8,9]+ abs.(M₁_int)[1,2,ℓ]*c₁^(ℓ-1)
    end
end
for ℓ in axes(M₂,3)
    if size(M₂,1) >= 2 
        if ℓ >= 2
            ΔDF₂₅[5,1] = ΔDF₂₅[5,1] + abs.(M₂_int)[2,1,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₂ 
            ΔDF₂₅[10,1] = ΔDF₂₅[10,1] + abs.(M₂_int)[2,1,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₉   
        end
        ΔDF₂₅[5,2] = ΔDF₂₅[5,2]+ abs.(M₂_int)[2,1,ℓ]*c₁^(ℓ-1)
        ΔDF₂₅[10,7] = ΔDF₂₅[5,2]+ abs.(M₂_int)[2,1,ℓ]*c₁^(ℓ-1)
    end
    if size(M₂,2) >= 2 && ℓ >= 2
        if ℓ >= 2
            ΔDF₂₅[5,1] = ΔDF₂₅[5,1] + abs.(M₂_int)[1,2,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₄  
            ΔDF₂₅[10,1] = ΔDF₂₅[10,1] + abs.(M₂_int)[1,2,ℓ]*interval(ℓ-1)*c₁^(ℓ-2)*c₂₅_₉  
        end
        ΔDF₂₅[5,4] = ΔDF₂₅[5,4]+ abs.(M₂_int)[1,2,ℓ]*c₁^(ℓ-1)
        ΔDF₂₅[10,9] = ΔDF₂₅[10,9]+ abs.(M₂_int)[1,2,ℓ]*c₁^(ℓ-1)
    end
end
ΔDF₂₅[3,1][0] = ExactReal(0)
ΔDF₂₅[3,2][0] = ExactReal(0)
ΔDF₂₅[3,4][0] = ExactReal(0)
ΔDF₂₅[5,1][0] = ExactReal(0)
ΔDF₂₅[5,2][0] = ExactReal(0)
ΔDF₂₅[5,4][0] = ExactReal(0)
ΔDF₂₅[8,1][0] = ExactReal(0)
ΔDF₂₅[8,7][0] = ExactReal(0)
ΔDF₂₅[8,9][0] = ExactReal(0)
ΔDF₂₅[10,1][0] = ExactReal(0)
ΔDF₂₅[10,7][0] = ExactReal(0)
ΔDF₂₅[10,9][0] = ExactReal(0)

ΔDF₂₅_temp = sum( ΔDF₂₅[:,2:end] , dims = 2)
ΔDF₂₅_₁ =  ΔDF₂₅[:,1] 

counter_order_poly_1 = [0;0]
counter_order_poly_2 = [0;0]
for i = 1:10
    local length_ΔDF₂₅_₁ = length(ΔDF₂₅_₁[i][:]) 
    local length_ΔDF₂₅_temp = length(ΔDF₂₅_temp[i][:]) 
    if length_ΔDF₂₅_temp > counter_order_poly_1[1]
        counter_order_poly_1[:] = [ length_ΔDF₂₅_temp , i]
    end
    if length_ΔDF₂₅_₁ > counter_order_poly_2[1]
        counter_order_poly_2[:] = [ length_ΔDF₂₅_₁ , i]
    end
end


ΔDF₂₅ = copy( ΔDF₂₅_temp[ counter_order_poly_1[2]  ] )
ΔDF₂₅_₁_ = copy( ΔDF₂₅_₁[ counter_order_poly_2[2]  ] )

for i = 1:10
    local length_ΔDF₂₅_temp = length(ΔDF₂₅_temp[i][:]) 
    if length_ΔDF₂₅_temp > 0
        for j = 1:length_ΔDF₂₅_temp
            ΔDF₂₅[j] = maximum( [ΔDF₂₅_temp[i][j] ; ΔDF₂₅[j]] )
        end
    end
    local length_ΔDF₂₅_₁ = length(ΔDF₂₅_₁[i][:]) 
    if length_ΔDF₂₅_₁ > 0
        for j = 1:length_ΔDF₂₅_₁
            ΔDF₂₅_₁_[j] = maximum( [ΔDF₂₅_₁[i][j] ; ΔDF₂₅_₁_[j]] )
        end
    end
end
ΔDF[25 , 25] = ΔDF₂₅
ΔDF[25 , 1] = ΔDF₂₅_₁_

vect_p = interval.([0;0;0;0]) 
vect_pₙ = interval.([0;0;0;0]) 
vect_pₘ = interval.([0;0;0;0]) 
for n = 1:4
    for i = 0:order(p_int)[n][1]
        for j = 0:order(p_int)[n][2]
            vect_p[n] = vect_p[n] +  interval(i) +  interval(j)
            vect_pₙ[n] = vect_pₙ[n] +  interval(i) 
            vect_pₘ[n] = vect_pₘ[n]  +  interval(j)
        end
    end
end
λ₁₁ = abs(eigenpairs_int[1])
λ₂₁ = abs(eigenpairs_int[6])

tails_p = interval(1)/minimum([ λ₁₁^2 , λ₁₁*λ₂₁, λ₂₁^2 ,  λ₂₁^2 ])
# i = 26
ΔDF[ 26 , 26 ] = Polynomial( [ interval(0) ;  ( vect_p[1] + tails_p ) / weight_x[25] ], :r)
ΔDF[ 26 , 25 ] = Polynomial( [ interval(0) ;  max(vect_pₙ, vect_pₘ ) / weight_x[26] ], :r)
# i = 27
ΔDF[ 27 , 27 ] = Polynomial( [ interval(0) ;  ( vect_p[2] + tails_p ) / weight_x[26] ], :r)
ΔDF[ 27 , 25 ] = Polynomial( [ interval(0) ;  max(vect_pₙ, vect_pₘ ) / weight_x[27] ], :r)
# i = 28
ΔDF[ 28 , 28 ] = Polynomial( [ interval(0) ;  ( vect_p[3] + tails_p ) / weight_x[25] ], :r)
ΔDF[ 28 , 25 ] = Polynomial( [ interval(0) ;  max(vect_pₙ, vect_pₘ ) / weight_x[28] ], :r)
# i = 29
ΔDF[ 29 , 29 ] = Polynomial( [ interval(0) ;  ( vect_p[4] + tails_p ) / weight_x[25] ], :r)
ΔDF[ 29 , 25 ] = Polynomial( [ interval(0) ;  max(vect_pₙ, vect_pₘ ) / weight_x[29] ], :r)

# Da loops
for n  in axes(M₁ ,1)
    for m in axes(M₁,2)
        for ℓ in axes(M₁,3)
            if M₁[n,m,ℓ] != 0
                ΔDF[ 3 , 24 ] = ΔDF[ 3 , 24 ] + (c₁^(ℓ-1))*abs.(M₁_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-1))
                ΔDF[ 12 , 24 ] = ΔDF[ 12 , 24 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ]*(c₁₁^(n-1)*c₁₃^(m-1))
                ΔDF[ 20 , 15 ] = ΔDF[ 20 , 15 ] + norm_T/interval(2)*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ]*(c₂₀^(n-1)*c₂₂^(m-1))
                if ℓ >= 2
                    ΔDF[ 3 , 1 ] = ΔDF[ 3 , 1 ] + interval((ℓ-1))*(c₁^(ℓ-2))*c₂₄*abs.(M₁_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-1))
                    ΔDF[ 12 , 1 ] = ΔDF[ 12 , 1 ] + interval(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₁_int)[n,m,ℓ] * c₁₁^(n-1)*c₁₃^(m-1)
                    ΔDF[ 27 , 1 ] = ΔDF[ 27 , 1 ] + interval(ℓ-1)*c₁^(ℓ-2)*abs.(M₁_int)[n,m,ℓ] * c₂₆^(n-1)*c₂₈^(m-1)
                    ΔDF[ 20 , 1 ] = ΔDF[ 20 , 1 ] + interval(ℓ-1)*c₁^(ℓ-2)*abs.(M₁_int)[n,m,ℓ] * c₂₀^(n-1)*c₂₂^(m-1)*norm_T/interval(2)*c₁₅
                end
                if n >= 2
                    ΔDF[ 3 , 2 ] = ΔDF[ 3 , 2 ] + interval((n-1))*(c₁^(ℓ-1))*c₂₄*abs.(M₁_int)[n,m,ℓ]*(c₂^(n-2)*c₄^(m-1))
                    ΔDF[ 12 , 11 ] = ΔDF[ 12 , 11 ] + c₁^(ℓ-1)*c₂₄*abs.(M₁_int)[n,m,ℓ]*interval(n-1)*(c₁₁^(n-2)*c₁₃^(m-1))
                    ΔDF[ 27 , 26 ] = ΔDF[ 27 , 26 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ]*interval(n-1)*(c₂₆^(n-2)*c₂₈^(m-1))
                    ΔDF[ 20 , 21 ] = ΔDF[ 20 , 21 ] + interval(n-1)*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * c₂₀^(n-2)*c₂₂^(m-1)*norm_T/interval(2)*c₁₅
                end
                if m >= 2
                    ΔDF[ 3 , 4 ] = ΔDF[ 3 , 4 ] + interval((m-1))*(c₁^(ℓ-1))*c₂₄*abs.(M₁_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-2))
                    ΔDF[ 12 , 13 ] = ΔDF[ 12 , 13 ] + c₁^(ℓ-1)*c₂₄*abs.(M₁_int)[n,m,ℓ]*interval(m-1)*(c₁₁^(n-1)*c₁₃^(m-2))
                    ΔDF[ 27 , 28 ] = ΔDF[ 27 , 28 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ]*interval(m-1)*(c₂₆^(n-1)*c₂₈^(m-2))
                    ΔDF[ 20 , 23 ] = ΔDF[ 20 , 23 ] + interval(m-1)*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * c₂₀^(n-1)*c₂₂^(m-2)*norm_T/interval(2)*c₁₅
                end
                if n-2 >= 0
                    ΔDF[ 8 , 7 ] = ΔDF[ 8 , 7 ]  + c₁^(ℓ-1)*c₂₄*abs.(M₁_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1) 
                end
                if m-2 >= 0
                    ΔDF[ 8 , 9 ] = ΔDF[ 8 , 9 ]  + c₁^(ℓ-1)*c₂₄*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2) 
                end
                if n-2 >= 0 && m-2>= 0
                    ΔDF[ 8 , 24 ] = ΔDF[ 8 , 24 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * (interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ + interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉)
                elseif n-2 >= 0
                    ΔDF[ 8 , 24 ] = ΔDF[ 8 , 24 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇  
                elseif m-2 >= 0
                    ΔDF[ 8 , 24 ] = ΔDF[ 8 , 24 ] + c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₇  
                end
                if n-3 >= 0 
                    ΔDF[ 8 , 2 ] = ΔDF[ 8 , 2 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(n-1)*interval(n-2)*c₂_ext^(n-3)*c₄_ext^(m-1)*c₇ 
                end
                if n-2 >= 0 && m-2 >= 0
                    ΔDF[ 8 , 2 ] = ΔDF[ 8 , 2 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-2)*c₉ 
                    ΔDF[ 8 , 4 ] = ΔDF[ 8 , 4 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-2)*c₇ 
                end
                if m-3 >=0
                    ΔDF[ 8 , 4 ] = ΔDF[ 8 , 4 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*(m-2)*c₂_ext^(n-1)*c₄_ext^(m-3)*c₉ 
                end
                if n-2 >= 0 && m-2>= 0 && ℓ >= 2
                    ΔDF[ 8 , 1 ] = ΔDF[ 8 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₁_int)[n,m,ℓ] * (interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ + interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉ )
                elseif n-2 >= 0 && ℓ >= 2
                    ΔDF[ 8 , 1 ] = ΔDF[ 8 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₁_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ 
                elseif m-2 >= 0  && ℓ >= 2
                    ΔDF[ 8 , 1 ] = ΔDF[ 8 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₁_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉  
                end
            end
        end
    end
end
for n  in axes(M₂ ,1)
    for m in axes(M₂,2)
        for ℓ in axes(M₂,3)
            if M₂[n,m,ℓ] != 0
                ΔDF[ 5 , 24 ] = ΔDF[ 5 , 24 ] + (c₁^(ℓ-1))*abs.(M₂_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-1))
                ΔDF[ 14 , 24 ] = ΔDF[ 14 , 24 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ]*(c₁₁^(n-1)*c₁₃^(m-1))
                ΔDF[ 22 , 15 ] = ΔDF[ 22 , 15 ] + norm_T/interval(2)*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ]*(c₂₀^(n-1)*c₂₂^(m-1))
                if ℓ >= 2
                    ΔDF[ 5 , 1 ] = ΔDF[ 5 , 1 ] + interval((ℓ-1))*(c₁^(ℓ-2))*c₂₄*abs.(M₂_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-1))
                    ΔDF[ 14 , 1 ] = ΔDF[ 14 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₂_int)[n,m,ℓ] * c₁₁^(n-1)*c₁₃^(m-1)
                    ΔDF[ 29 , 1 ] = ΔDF[ 29 , 1 ] + interval(ℓ-1)*c₁^(ℓ-2)*abs.(M₂_int)[n,m,ℓ] * c₂₆^(n-1)*c₂₈^(m-1)
                    ΔDF[ 22 , 1 ] = ΔDF[ 22 , 1 ] + interval(ℓ-1)*c₁^(ℓ-2)*abs.(M₂_int)[n,m,ℓ] * c₂₀^(n-1)*c₂₂^(m-1)*norm_T/interval(2)*c₁₅
                end
                if n >= 2
                    ΔDF[ 5 , 2 ] = ΔDF[ 5 , 2 ] + interval((n-1))*(c₁^(ℓ-1))*c₂₄*abs.(M₂_int)[n,m,ℓ]*(c₂^(n-2)*c₄^(m-1))
                    ΔDF[ 14 , 11 ] = ΔDF[ 14 , 11 ] + c₁^(ℓ-1)*c₂₄*abs.(M₂_int)[n,m,ℓ]*interval(n-1)*(c₁₁^(n-2)*c₁₃^(m-1))
                    ΔDF[ 29 , 26 ] = ΔDF[ 29 , 26 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ]*interval(n-1)*(c₂₆^(n-2)*c₂₈^(m-1))
                    ΔDF[ 22 , 21 ] = ΔDF[ 22 , 21 ] + interval(n-1)*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * c₂₀^(n-2)*c₂₂^(m-1)*norm_T/interval(2)*c₁₅
                end
                if m >= 2
                    ΔDF[ 5 , 4 ] = ΔDF[ 5 , 4 ] + interval((m-1))*(c₁^(ℓ-1))*c₂₄*abs.(M₂_int)[n,m,ℓ]*(c₂^(n-1)*c₄^(m-2))
                    ΔDF[ 14 , 13 ] = ΔDF[ 12 , 13 ] + c₁^(ℓ-1)*c₂₄*abs.(M₂_int)[n,m,ℓ]*interval(m-1)*(c₁₁^(n-1)*c₁₃^(m-2))
                    ΔDF[ 29 , 28 ] = ΔDF[ 29 , 28 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ]*interval(m-1)*(c₂₆^(n-1)*c₂₈^(m-2))
                    ΔDF[ 22 , 23 ] = ΔDF[ 22 , 23 ] + interval(m-1)*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * c₂₀^(n-1)*c₂₂^(m-2)*norm_T/interval(2)*c₁₅
                end
                if n-2 >= 0
                    ΔDF[ 10 , 7 ] = ΔDF[ 10 , 7 ]  + c₁^(ℓ-1)*c₂₄*abs.(M₂_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1) 
                end
                if m-2 >= 0
                    ΔDF[ 10 , 9 ] = ΔDF[ 10 , 9 ]  + c₁^(ℓ-1)*c₂₄*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2) 
                end
                if n-2 >= 0 && m-2>= 0
                    ΔDF[ 10 , 24 ] = ΔDF[ 10 , 24 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * (interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ + (m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉)
                elseif n-2 >= 0
                    ΔDF[ 10 , 24 ] = ΔDF[ 10 , 24 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇
                elseif m-2 >= 0
                    ΔDF[ 10 , 24 ] = ΔDF[ 10 , 24 ] + c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₇
                end
                if n-3 >= 0 
                    ΔDF[ 10 , 2 ] = ΔDF[ 10 , 2 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(n-1)*interval(n-2)*c₂_ext^(n-3)*c₄_ext^(m-1)*c₇ 
                end
                if n-2 >= 0 && m-2 >= 0
                    ΔDF[ 10 , 2 ] = ΔDF[ 10 , 2 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-2)*c₉ 
                    ΔDF[ 10 , 4 ] = ΔDF[ 10 , 4 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-2)*c₇ 
                end
                if m-3 >=0
                    ΔDF[ 10 , 4 ] = ΔDF[ 10 , 4 ]  + c₂₄*c₁^(ℓ-1)*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*interval(m-2)*c₂_ext^(n-1)*c₄_ext^(m-3)*c₉ 
                end
                if n-2 >= 0 && m-2>= 0 && ℓ >= 2
                    ΔDF[ 10 , 1 ] = ΔDF[ 10 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₂_int)[n,m,ℓ] * (interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ + interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉ )
                elseif n-2 >= 0 && ℓ >= 2
                    ΔDF[ 10 , 1 ] = ΔDF[ 10 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₂_int)[n,m,ℓ] * interval(n-1)*c₂_ext^(n-2)*c₄_ext^(m-1)*c₇ 
                elseif m-2 >= 0  && ℓ >= 2
                    ΔDF[ 10 , 1 ] = ΔDF[ 10 , 1 ] + ExactReal(ℓ-1)*c₁^(ℓ-2)*c₂₄*abs.(M₂_int)[n,m,ℓ] * interval(m-1)*c₂_ext^(n-1)*c₄_ext^(m-2)*c₉ 
                end
            end
        end
    end
end
w_vect = Vector{Any}(undef , 4)
for n = 1:4
    w_vect[n] =  zeros(eltype(w_int), space(w)[n])
    for i = -order(w)[n][1]:order(w)[n][1]
        for j = 2:order(w)[n][2]
            w_vect[n][(i,j)] = interval(j) / ν_w^(abs(interval(i)) + interval(j))
        end
    end
end
B = Matrix{Any}(undef,29,29)
for i = 1:29
    for j = 1:29
        for k = 1:29
            if isassigned(B,i,j)
                B[i,j] = B[i,j] + A_[i,k]*ΔDF[k,j]
            else
                B[i,j] =  A_[i,k]*ΔDF[k,j]
            end
        end
    end
end

for i = 1:29
    if i in [1;6;15;16;17;18;19;24;25]
        local space_test_in = ℓ∞()
    elseif i in [ 2 ;3;4;5]
        local space_test_in =  ℓ¹(GeometricWeight(νᵧ))
    elseif i in [ 7;8;9;10]
        local space_test_in =  ℓ¹(GeometricWeight(νᵥ))
    elseif i in [ 11;12;13;14]
        local space_test_in =  ℓ¹(GeometricWeight(ν_w))
    elseif i in [ 20;21;22;23]
        local space_test_in =  ℓ¹(GeometricWeight(νₐ))
    elseif i in [ 26;27;28;29]
        local space_test_in =  ℓ¹(GeometricWeight(νₚ))
    end
    for n = 1:4
        temp_B = norm( Sequence( space_X[i] ,  norm.( Aᵢⱼ_cheb[i,10+n] , 1 )*w_vect[n][:]) , space_test_in )
        if i == 10+n
            temp_B = max( temp_B , norm( interval(1) / λ_cheb[1](interval( -1,1 )) ,1) ) 
        end
        B[i,6] = B[i,6] + Polynomial( [ interval(0) , temp_B/weight_x[10+n] ]  ,:r) 
        B[i,10+n] = B[i,10+n] + Polynomial( [ interval(0) , temp_B/weight_x[6] ]  ,:r) 
    end
    Aᵢ₁₉_cheb = norm( Sequence( space_X[i] , norm.( Aᵢⱼ_cheb[i,19] , 1 )[:,1] ) , space_test_in ) 
    Aᵢ₂₀_cheb = norm( Sequence( space_X[i] , norm.( Aᵢⱼ_cheb[i,20] , 1 )[:,1] ) , space_test_in ) 
    Aᵢ₂₁_cheb = norm( Sequence( space_X[i] , norm.( Aᵢⱼ_cheb[i,21] , 1 )[:,1] ) , space_test_in ) 
    Aᵢ₂₂_cheb = norm( Sequence( space_X[i] , norm.( Aᵢⱼ_cheb[i,22] , 1 )[:,1] ) , space_test_in ) 
    B[i,16] = B[i,17] + Aᵢ₁₉_cheb *ΔDF_₁₉_₁₆ + Aᵢ₂₀_cheb *ΔDF_₂₀_₁₆ + Aᵢ₂₁_cheb *ΔDF_₂₁_₁₆ + Aᵢ₂₂_cheb *ΔDF_₂₂_₁₆
    B[i,17] = B[i,17] + Aᵢ₁₉_cheb *ΔDF_₁₉_₁₇ + Aᵢ₂₀_cheb *ΔDF_₂₀_₁₇ + Aᵢ₂₁_cheb *ΔDF_₂₁_₁₇ + Aᵢ₂₂_cheb *ΔDF_₂₂_₁₇ 
end

Z₂r_temp = sum( B .* (weight_x*transpose( interval(1) ./ weight_x )) ,dims = 2)
for i =1:29
    Z₂r_temp[i][0] =ExactReal( 0 )
end

counter_order_poly = [0;0]
for i = 1:29
    local length_Z₂rⁱ = length(Z₂r_temp[i][:]) 
    if length_Z₂rⁱ > counter_order_poly[1]
        counter_order_poly[:] = [ length_Z₂rⁱ , i]
    end
end

Z₂r = copy( Z₂r_temp[ counter_order_poly[2]  ] )
for i = 1:29
    local length_Z₂rⁱ = length(Z₂r_temp[i][:]) 
    if length_Z₂rⁱ > 0
        for j = 1:length_Z₂rⁱ
            Z₂r[j] = maximum( [Z₂r_temp[i][j] ; Z₂r[j]] )
        end
    end
end
printstyled("✶ Z₂ = $(Polynomial( sup.(Z₂r[1:end]) ,:r ))\n"; color = :green)
# =========================================================================
#  Radii Polynomial
# ========================================================================= 
printstyled("============================================================\n"; color = :blue)
printstyled("5. Roots of Radii Polynomial\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)

poly_int =  Z₂r*Polynomial( [interval(0), interval(1)] ,:r ) + Polynomial( [ Y₀ , -( 1 - Z₁ - Z₀) ] ,:r )

poly = Polynomial( sup.(poly_int[:]) )
rootsy = roots(poly)
root_all = [0.0;0.0]
counter = 1
r_star = min( (νₚ - norm_z)/(norm_t*interval(10))  , (ν_w - norm_y)/(norm_h*interval(10))  )
for i in axes(rootsy,1)
    if imag(rootsy[i]) == 0 && real(rootsy[i]) > 0
        root_all[counter] =  real(rootsy[i])
        counter = counter + 1
    end
end
if isempty(root_all)
    error("Proof did not work no real positive roots")
elseif root_all[1] > inf(r_star)
    error(string("Smallest positive root is bigger than r_star =$(inf(r_star) ) "))
else
    printstyled("✶ r_min = $(root_all[1])\n"; color = :green)
    printstyled("✶ r_max = $(root_all[2])\n"; color = :green)
end
Y₀ = sup(Y₀);
Z₀ = sup(Z₀);
Z₁ = sup(Z₁);
printstyled("============================================================\n"; color = :blue)
printstyled("Proof Complete\n"; color = :blue)
printstyled("============================================================\n"; color = :blue)

    return X_cheb , η_cheb, X_grid , η_grid, N, κ , ξ₁ᵏ, ξ₂ᵏ, star, lengthᵥ, M₁, M₂ , Y₀ , Z₀, Z₁ , Z₂r , root_all , νᵧ, νᵥ , νₐ , ν_w ,  νₚ , weight_x
end
#  dir_name = joinpath(@__DIR__, "Data", "GS Orientable", "PROVEN","Data_GS_snaking_proof_")
# jldsave(dir_name*"14.jld2"; X_cheb , η_cheb, X_grid , η_grid, N, κ , ξ₁ᵏ, ξ₂ᵏ, star, lengthᵥ, M₁, M₂ , Y₀ , Z₀, Z₁ , Z₂r , root_all , νᵧ, νᵥ , νₐ , ν_w ,  νₚ , weight_x)
# Need to redo 1-7


