# ==============================================================
# Finite difference approximation of derivative of the maps
# ==============================================================

function finite_diff_γ_Γ!(DF,γ,ω,M₁,M₂,α) 
    N = size(DF,2)
    E =  1*convert.(eltype(γ),Matrix(I,size(DF,1),N))
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(γ),  E[1:end,i]   ) 
        F₊ = Γ!(Sequence( space(Derivative(1)*γ)  , zeros(eltype(γ),length(γ),1)[:] ),γ+e,ω,M₁,M₂,α) 
        F₋ = Γ!(Sequence( space(Derivative(1)*γ)  , zeros(eltype(γ),length(γ),1)[:] ),γ-e,ω,M₁,M₂,α) 
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end
function finite_diff_ω_Γ!(DF,γ,ω,M₁,M₂,α) 
    e = 10^(-6)
    F₊ = Γ!(zeros( eltype(γ) , space(Derivative(1)*γ) ) ,γ,ω+e,M₁,M₂,α) 
    F₋ = Γ!(zeros( eltype(γ) ,space(Derivative(1)*γ) ) ,γ,ω-e,M₁,M₂,α) 
    DF[:] = (F₊[:] - F₋[:])/(2*e)

    return DF
end
function finite_diff_V!( DF, ω, v, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ)
    N = length(v)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 0:N
        if i == 0
            e = 10^(-6)
            F₊ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v, γ_ext, M₁, M₂,α, λ+e, κ, lengthᵥ)
            F₋ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v, γ_ext, M₁, M₂,α, λ-e, κ, lengthᵥ)
            DF[:,1] = (F₊[:] - F₋[:])/(2*h)
        else
            e = h*Sequence(space(v),  E[1:end,i]   ) 
            F₊ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v+e, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ)
            F₋ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v-e, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ) 
            DF[:,i+1] = (F₊[:] - F₋[:])/(2*h)
        end

    end
    return DF
end

function finite_diff_ω_V!( DF, ω, v, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ)
    e = 10^(-6)
    F₊ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω+e, v, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ)
    F₋ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω-e, v, γ_ext, M₁, M₂,α, λ, κ, lengthᵥ)
    DF[:] = (F₊[:] - F₋[:])/(2*e)
    return DF
end
function finite_diff_E!(DF,eigenpairs, M₁, M₂, ξ₁ᵏ ,ξ₂ᵏ, star) 
    N = length(eigenpairs)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(eigenpairs),  E[1:end,i]   ) 
        F₊ = E!(zero(eigenpairs),eigenpairs+e, M₁, M₂, ξ₁ᵏ ,ξ₂ᵏ, star) 
        F₋ = E!(zero(eigenpairs),eigenpairs-e, M₁, M₂, ξ₁ᵏ ,ξ₂ᵏ, star)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_ω_W!(DF, w, γ_ext, v, ω, λ, M₁, M₂,α)
    e = 10^(-6)
    F₊ = W!(zero(w), w, γ_ext, v, ω+e, λ, M₁, M₂,α)
    F₋ = W!(zero(w), w, γ_ext, v, ω-e, λ, M₁, M₂,α)
    DF[:] = (F₊[:] - F₋[:])/(2*e)
    return DF
end

function finite_diff_λ_W!(DF, w, γ_ext, v, ω, λ, M₁, M₂,α)
    e = 10^(-6)
    F₊ = W!(zero(w), w, γ_ext, v, ω, λ+e, M₁, M₂,α)
    F₋ = W!(zero(w), w, γ_ext, v, ω, λ-e, M₁, M₂,α)
    DF[:] = (F₊[:] - F₋[:])/(2*e)
    return DF
end

function finite_diff_v_W!(DF, w, γ_ext, v, ω, λ, M₁, M₂,α)
    N = length(v)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(v),  E[1:end,i]   ) 
        F₊ = W!(zero(w), w, γ_ext, v+e, ω, λ, M₁, M₂,α)
        F₋ = W!(zero(w), w, γ_ext, v-e, ω, λ, M₁, M₂,α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_w_W!(DF, w, γ_ext, v, ω, λ, M₁, M₂,α)
    N = length(w)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(w),  E[1:end,i]   ) 
        F₊ = W!(zero(w), w+e, γ_ext, v, ω, λ, M₁, M₂,α)
        F₋ = W!(zero(w), w-e, γ_ext, v, ω, λ, M₁, M₂,α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end
function finite_diff_e_P!(DF,p, eigenpairs, M₁, M₂,α) 
    N = length(eigenpairs)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(eigenpairs),  E[1:end,i]   ) 
        F₊ = P!(zero(p), p, eigenpairs+e, M₁, M₂,α)
        F₋ = P!(zero(p), p, eigenpairs-e, M₁, M₂,α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end
function finite_diff_p_P!(DF,p, eigenpairs, M₁, M₂,α) 
    N = length(p)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:length(project(zero(p), domain(DF) ))
        e = h*Sequence(space(p),  E[1:end,i]   ) 
        F₊ = P!(zero(p) , p+e, eigenpairs, M₁, M₂,α)
        F₋ = P!(zero(p), p-e, eigenpairs, M₁, M₂,α)
        DF[:,i] = project((F₊[:] - F₋[:])/(2*h),space(p))
    end
    return DF
end

function finite_diff_a_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    N = length(a)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()× ParameterSpace()) 
    for i = 1:N
        e = h*Sequence(space(a),  E[1:end,i]   ) 
        F₊ = G!(copy(G), a+e, w, p, θ₀, L, σ, M₁, M₂, α)
        F₋ = G!(copy(G), a-e, w, p, θ₀, L, σ, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_θ₀_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    N = length(σ)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()× ParameterSpace()) 
    for i = 1:N
        e = h*Sequence(space(σ),  E[1:end,i]   ) 
        F₊ = G!(copy(G), a, w, p, θ₀+e, L, σ, M₁, M₂, α)
        F₋ = G!(copy(G), a, w, p, θ₀-e, L, σ, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_L_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()) 
    e = h
    F₊ = G!(copy(G), a, w, p, θ₀, L+e, σ, M₁, M₂, α)
    F₋ = G!(copy(G), a, w, p, θ₀, L-e, σ, M₁, M₂, α)
    DF[:] = (F₊[:] - F₋[:])/(2*h)

    return DF
end

function finite_diff_σ_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    N = length(σ)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()) 
    for i = 1:N
        e = h*Sequence(space(σ),  E[1:end,i]   ) 
        F₊ = G!(copy(G), a, w, p, θ₀, L, σ+e, M₁, M₂, α)
        F₋ = G!(copy(G), a, w, p, θ₀, L, σ-e, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_p_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    N = length(p)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()) 
    for i = 1:N
        e = h*Sequence(space(p),  E[1:end,i]   ) 
        F₊ = G!(copy(G), a, w, p+e, θ₀, L, σ, M₁, M₂, α)
        F₋ = G!(copy(G), a, w, p-e, θ₀, L, σ, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_w_G!(DF, a, w, p, θ₀, L, σ, M₁, M₂, α)
    N = length(w)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    G = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()× ParameterSpace()) 
    for i = 1:N
        e = h*Sequence(space(w),  E[1:end,i]   ) 
        F₊ = G!(copy(G), a, w+e, p, θ₀, L, σ, M₁, M₂, α)
        F₋ = G!(copy(G), a, w-e, p, θ₀, L, σ, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_γ_V!( DF, ω, v, γ, M₁, M₂, α, λ, κ, lengthᵥ)
    N = length(γ)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(γ),  E[1:end,i]   ) 
        F₊ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v, γ+e, M₁, M₂, α, λ, κ, lengthᵥ)
        F₋ = V!( Sequence( ParameterSpace() × space(v),zeros(eltype(v),length(v)+1,1)[:] ), ω, v, γ-e, M₁, M₂, α, λ, κ, lengthᵥ)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end
function finite_diff_α_F!(DF, x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
    h = 10^(-6)
    F = zeros(ComplexF64,space(Derivative(1)*component(x,1:4))×space(x)[5:27])
    e = h
    F₊ = F!(copy(F), x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α+e)
    F₋ = F!(copy(F), x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α-e)
    DF[:] = (F₊[:] - F₋[:])/(2*h)

    return DF
end
function finite_diff_x_F!(DF, x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
    N = length(x)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    F = zeros(ComplexF64,space_f!(x))
    for i = 1:N
        e = h*Sequence(space(x),  E[1:end,i]   ) 
        F₊ = F!(copy(F), x+e, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
        F₋ = F!(copy(F), x-e, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end

function finite_diff_γ_W!(DF, w, γ_ext, v, ω, λ, M₁, M₂,α)
    N = length(γ)
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = h*Sequence(space(γ),  E[1:end,i]   ) 
        F₊ = W!(zero(w), w, γ_ext+e, v, ω, λ, M₁, M₂,α)
        F₋ = W!(zero(w), w, γ_ext-e, v, ω, λ, M₁, M₂,α)
        DF[:,i] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end


function Dₐ₁G₅_fin!Fᵢ!(DF, X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , η  )
    α = real(X[1])
    x = x_remove_complex!(component(X,2:29))
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
    a₁,a₂,a₃,a₄ = eachcomponent(a)

    N = length(zeros(domain(DF)))
    E =  1*Matrix(I,N,N)
    h = 10^(-6)
    for i = 1:N
        e = Sequence( domain(DF)  ×  space(a₂) × space(a₃) × space(a₄) ,  [ h*E[1:end,i]  ; zero(a₂)[:] ; zero(a₃)[:] ; zero(a₄)[:] ])

        F₊ = Gᵢ!( zeros(eltype(p),codomain(DF)) , a+e , w, p, θ₀, L, σ, M₁, M₂, α , 8 )
        F₋ = Gᵢ!( zeros(eltype(p),codomain(DF)) , a-e , w, p, θ₀, L, σ, M₁, M₂, α , 8 )
        DF[:,i-1] = (F₊[:] - F₋[:])/(2*h)
    end
    return DF
end


