# ====================================================
# Maps
# ====================================================

# Periodic Solution
function Γ!(Γ,γ,ω,M₁,M₂,α) 
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
    D = Derivative(1)

    Γ₁ = project(-(D*γ₁) + ω*γ₂ , space(Γ)[1])
    Γ₂ = project((D*γ₂), space(Γ)[2])
    Γ₃ = project(-(D*γ₃) + ω*γ₄, space(Γ)[3])
    Γ₄ = project((D*γ₄), space(Γ)[4])

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    Γ₂ = Γ₂ - project(α^(ℓ-1)*ω*M₁[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(Γ)[2])
                end
            end
        end
    end
    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    Γ₄ = Γ₄ - project(α^(ℓ-1)*ω*M₂[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(Γ)[4])
                end
            end
        end
    end

    Γ[:] = [ Γ₁[:] ; Γ₂[:] ; Γ₃[:] ; Γ₄[:] ]
    return Γ
end
# Bundle of periodic solution
function V!( V, ω, v, γ, M₁, M₂,α, λ, κ, lengthᵥ)
    γ_ext =  γ_extened_orientable!(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    v_κ = component(v,κ)

    γ₁ = component(γ_ext,1)
    γ₃ = component(γ_ext,3)

    D = Derivative(1)
    
    V₀ =  evaluate(v_κ,0) - lengthᵥ
    V₁ = project((D*v₁) + λ*v₁ - ω*v₂,space(V)[2])
    V₂ = project((D*v₂) + λ*v₂,space(V)[3])
    V₃ = project((D*v₃) + λ*v₃ - ω*v₄,space(V)[4])
    V₄ = project((D*v₄) + λ*v₄ ,space(V)[5])

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n-2 >= 0 && m-2>= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[3])
                    elseif n-2 >= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ ,space(V)[3])
                    elseif m-2 >= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[3])
                    end
                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                 if M₂[n,m,ℓ] != 0
                    if n-2 >= 0 && m-2>= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃  ,space(V)[5])
                    elseif n-2 >= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁  ,space(V)[5])
                    elseif m-2 >= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[5])
                    end
                end
            end
        end
    end

    V[:] = [V₀; V₁[:] ; V₂[:] ; V₃[:] ; V₄[:]] 
    return V
end
# Local Unstable manifold of periodic solution
function W!(W, w, γ, v, ω, λ, M₁, M₂,α)
    γ_ext =  γ_extened_orientable!(γ)
    w₁,w₂,w₃,w₄ = eachcomponent(w)
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ_ext)
    v₁,v₂,v₃,v₄ = eachcomponent(v)

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])

    W₁ = project(D₁₀*w₁ + λ*M*(D₀₁*w₁) - ω*w₂, space(component(W,1)))
    W₂ = project(D₁₀*w₂ + λ*M*(D₀₁*w₂) , space(component(W,2)))
    W₃ = project(D₁₀*w₃ + λ*M*(D₀₁*w₃) - ω*w₄, space(component(W,3)))
    W₄ = project(D₁₀*w₄ + λ*M*(D₀₁*w₄) , space(component(W,4)))

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    W₂ = W₂ - α^(ℓ-1)*ω*M₁[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(component(W,2)))
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    W₄ = W₄ - α^(ℓ-1)*ω*M₂[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(component(W,4)))
                end
            end
        end
    end

    W₁[(:,0)] = project(Sequence( space(w₁)[1], w₁[(:,0)] ) - γ₁,space(W₁)[1])[:]
    W₁[(:,1)] = project(Sequence( space(w₁)[1], w₁[(:,1)] ) - v₁,space(W₁)[1])[:]
    W₂[(:,0)] = project(Sequence( space(w₂)[1], w₂[(:,0)] ) - γ₂,space(W₂)[1])[:]
    W₂[(:,1)] = project(Sequence( space(w₂)[1], w₂[(:,1)] ) - v₂,space(W₂)[1])[:]
    W₃[(:,0)] = project(Sequence( space(w₃)[1], w₃[(:,0)] ) - γ₃,space(W₃)[1])[:]
    W₃[(:,1)] = project(Sequence( space(w₃)[1], w₃[(:,1)] ) - v₃,space(W₃)[1])[:]
    W₄[(:,0)] = project(Sequence( space(w₄)[1], w₄[(:,0)] ) - γ₄,space(W₄)[1])[:]
    W₄[(:,1)] = project(Sequence( space(w₄)[1], w₄[(:,1)] ) - v₄,space(W₄)[1])[:]

    W[:] =  [W₁[:]; W₂[:]; W₃[:]; W₄[:]] 
    return W
end
# Eigenpair of Equilibrium at the origin
function E!(E,eigenpairs, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star)
    if length(star) == 1
        star = [star; star]
    end
    λ₁ = component(eigenpairs,1)[1]
    ξ₁ = component(eigenpairs,2:5)[:]
    λ₂ = component(eigenpairs,6)[1]
    ξ₂ = component(eigenpairs,7:10)[:]

    if size(M₁,1) < 2 || size(M₁,2) < 2
        M₁_temp = zeros(typeof(M₁[1]),2,2,size(M₁,3))
        if size(M₁,1) >= 2 
            M₁_temp[2,1,:] = M₁[2,1,:] 
        end
        if size(M₁,2) >= 2 
            M₁_temp[1,2,:] = M₁[1,2,:] 
        end
    else
        M₁_temp = M₁
    end

    if size(M₂,1) < 2 || size(M₂,2) < 2
        M₂_temp = zeros(typeof(M₂[1]),2,2,size(M₂,3))
        if size(M₂,1) >= 2 
            M₂_temp[2,1,:] = M₂[2,1,:] 
        end
        if size(M₂,2) >= 2 
            M₂_temp[1,2,:] = M₂[1,2,:] 
        end
    else
        M₂_temp = M₂
    end
    if eltype(α) == Interval{Float64}
        M = interval.( [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0])
    else
        M = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    end
    for  ℓ in axes(M₁,3)
        M = M + [ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); α^ExactReal(ℓ-1)*M₁_temp[2,1,ℓ] ExactReal(0) α^ExactReal(ℓ-1)*M₁_temp[1,2,ℓ] ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0)]
    end
    for  ℓ in axes(M₂,3)
        M = M + [ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); α^ExactReal(ℓ-1)*M₂_temp[2,1,ℓ] ExactReal(0) α^ExactReal(ℓ-1)*M₂_temp[1,2,ℓ] ExactReal(0)]
    end

    E[1] = ξ₁[star[1]] - ξ₁ᵏ
    E[6]  = ξ₂[star[2]] - ξ₂ᵏ

    E[2:5] = λ₁*ξ₁ - M*ξ₁
    E[7:10] = λ₂*ξ₂ - M*ξ₂
    return E
end
# Local Stable Manifold of Equilibrium at the origin
function P!(P, p, eigenpairs, M₁, M₂, α)
    p₁, p₂, p₃, p₄ = eachcomponent(p)
    λ₁ = component(eigenpairs,1)[1]
    ξ₁ = component(eigenpairs,2:5)[:]
    λ₂ = component(eigenpairs,6)[1]
    ξ₂ = component(eigenpairs,7:10)[:] 

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
    M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])

    P₁ = project(λ₁*M₁₀*(D₁₀*p₁) + λ₂*M₀₁*(D₀₁*p₁) - p₂ , space(P)[1])
    P₂ = project(λ₁*M₁₀*(D₁₀*p₂) + λ₂*M₀₁*(D₀₁*p₂) , space(P)[2])
    P₃ = project(λ₁*M₁₀*(D₁₀*p₃) + λ₂*M₀₁*(D₀₁*p₃) - p₄ , space(P)[3])
    P₄ = project(λ₁*M₁₀*(D₁₀*p₄) + λ₂*M₀₁*(D₀₁*p₄) , space(P)[4])

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    P₂ = P₂ - M₁[i,j,ℓ] * α^(ℓ-1)*project(p₁^(i-1)*p₃^(j-1), space(P)[2])
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    P₄ = P₄ - M₂[i,j,ℓ] * α^(ℓ-1)*project(p₁^(i-1)*p₃^(j-1), space(P)[4])
                end
            end
        end
    end

    P₁[(0,0)] = p₁[(0,0)] 
    P₁[(1,0)] = p₁[(1,0)] - ξ₁[1]
    P₁[(0,1)] = p₁[(0,1)] - ξ₂[1]

    P₂[(0,0)] = p₂[(0,0)] 
    P₂[(1,0)] = p₂[(1,0)] - ξ₁[2]
    P₂[(0,1)] = p₂[(0,1)] - ξ₂[2]

    P₃[(0,0)] = p₃[(0,0)] 
    P₃[(1,0)] = p₃[(1,0)] - ξ₁[3]
    P₃[(0,1)] = p₃[(0,1)] - ξ₂[3]

    P₄[(0,0)] = p₄[(0,0)] 
    P₄[(1,0)] = p₄[(1,0)] - ξ₁[4]
    P₄[(0,1)] = p₄[(0,1)] - ξ₂[4]

    P[:] =  [P₁[:]; P₂[:]; P₃[:]; P₄[:]] 
    return P
end
# Connecting Orbit
function G!(G, a, w, p, θ, L, σ, M₁, M₂, α)
    a₁, a₂, a₃, a₄ = eachcomponent(a) 
    p₁, p₂, p₃, p₄ = eachcomponent(p) 
    w₁, w₂, w₃, w₄ = eachcomponent(w) 
    σ₁, σ₂ = eachcomponent(σ)
    θ₁, θ₂ = eachcomponent(θ)

    component(G,1)[1] = p₁(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₁(1)
    component(G,2)[1] = p₂(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₂(1)
    component(G,3)[1] = p₃(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₃(1)
    component(G,4)[1] = p₄(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₄(1)

    component(G,5)[:] = project( 2*(0:order(a₁)).*a₁ , space(G)[5] )[:]
    component(G,6)[:] = project( 2*(0:order(a₂)).*a₂ , space(G)[6] )[:]
    component(G,7)[:] = project( 2*(0:order(a₃)).*a₃ , space(G)[7] )[:]
    component(G,8)[:] = project( 2*(0:order(a₄)).*a₄ , space(G)[8] )[:]
  
    Ψ₁ = L/2*project( a₂ ,Chebyshev(order(component(G,5))+1) )
    Ψ₂ = Sequence(Chebyshev(order(component(G,6))+1),zeros(order(component(G,6))+2))
    Ψ₃ = L/2*project( a₄ ,Chebyshev(order(component(G,7))+1) )
    Ψ₄ = Sequence(Chebyshev(order(component(G,8))+1),zeros(order(component(G,8))+2))

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ  in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    Ψ₂ = Ψ₂ + L/2*M₁[i,j,ℓ] * α^(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₂))
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ  in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    Ψ₄ = Ψ₄ + L/2*M₂[i,j,ℓ] * α^(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₄))
                end
            end 
        end
    end

    component(G,5)[1:end] = component(G,5)[1:end] + Ψ₁[2:end]  - Ψ₁[0:end-2] 
    component(G,6)[1:end] = component(G,6)[1:end] + Ψ₂[2:end]  - Ψ₂[0:end-2] 
    component(G,7)[1:end] = component(G,7)[1:end] + Ψ₃[2:end]  - Ψ₃[0:end-2] 
    component(G,8)[1:end] = component(G,8)[1:end] + Ψ₄[2:end]  - Ψ₄[0:end-2] 

    N_w = order(w₁)    
    n = Sequence(space(w₁),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    m = Sequence(space(w₁),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    component(G,5)[0] = sum(w₁[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₁[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₁[(0,:)].*(-0.95).^m[(0,:)])
    component(G,5)[0] = component(G,5)[0]  - a₁(-1)

    N_w = order(w₂)    
    n = Sequence(space(w₂),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    m = Sequence(space(w₂),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    component(G,6)[0] = sum(w₂[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₂[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₂[(0,:)].*(-0.95).^m[(0,:)])
    component(G,6)[0] =  component(G,6)[0] - a₂(-1)

    N_w = order(w₃)    
    n = Sequence(space(w₃),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    m = Sequence(space(w₃),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    component(G,7)[0] = sum(w₃[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₃[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₃[(0,:)].*(-0.95).^m[(0,:)])
    component(G,7)[0] = component(G,7)[0] - a₃(-1)

    N_w = order(w₄)    
    n = Sequence(space(w₄),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    m = Sequence(space(w₄),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
    component(G,8)[0] = sum(w₄[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₄[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₄[(0,:)].*(-0.95).^m[(0,:)])
    component(G,8)[0] = component(G,8)[0] - a₄(-1)

    component(G,9)[1] = 0.95 - σ₁[1]^2 - σ₂[1]^2
    component(G,10)[1] = 1 - θ₁[1]^2 - θ₂[1]^2
    return G
end
# Zero finding problem f
function F!(F, x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
    # Components
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)

    # Periodic Solutions Γ
    component(F,1:4)[:] = Γ!(real(component(F,1:4)),γ,ω,M₁,M₂,α)[:] 

    # Bundle of the periodic solution V
    component(F,5:9)[:] = V!(component(F,5:9), ω, v, γ, M₁, M₂, α, λ, κ, lengthᵥ)[:]

    # Local Unstable Manifold of the periodic solution W
    component(F,10:13)[:] = W!(component(F,10:13), w, γ, v, ω, λ, M₁, M₂, α)[:]

    # Connecting orbit G
    component(F,14:23)[:] = G!(component(F,14:23), a, w, p, θ₀, L, σ, M₁, M₂, α)[:]

    # Eigenpair of the stable eigenvalue E
    component(F,24)[:] = E!(component(F,24),eigenpairs, M₁, M₂, α, ξ₁ᵏ ,ξ₂ᵏ, star)[:]

    # Local Stable Manifold of the equilibrium at the origin P
    component(F,25:28)[:] = P!(component(F,25:28), p, eigenpairs, M₁, M₂, α)[:]

    return F

end
# Zero finding problem f
function F_all!(F, X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot )
    α = real(X[1])
    #x = x_remove_complex!(component(X,2:29))
    x = component(X,2:29)
    component(F,2:29)[:] = F!(component(F,2:29), x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)[:]
    #F[1] = dot((X[:]- Xₖ[:]),Xₖ_dot)
    F[1] = dot(Xₖ_dot,(X[:]- Xₖ[:]))
    return F

end


function Γ_all!(F, X, ω, M₁, M₂ ,Xₖ , Xₖ_dot )
    α = real(component(X,1)[1])
    γ = component(X,2:5)

    component(F,2:5)[:] = Γ!( component(F,2:5) , γ, ω , M₁ , M₂ , α )[:]
    F[1] = dot(Xₖ_dot,X[:]- Xₖ[:])
    return F

end

# Map F component-wise
function Fᵢ!(F, X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot , i )
    α = real(X[1])
    x = x_remove_complex!(component(X,2:29))
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)

    if i == 1
        F[1] = dot(Xₖ_dot,(X[:]- Xₖ[:]))[1]
    elseif i in [2;3;4;5]
        F[:] = Γᵢ!(F,γ,ω,M₁,M₂,α , i - 1  )[:] 
    elseif i in [6;7;8;9;10]
        F[:] = Vᵢ!(F, ω, v, γ, M₁, M₂, α, λ, κ, lengthᵥ , i - 5  )[:]
    elseif i in [11;12;13;14]
        F[:] = Wᵢ!(F, w, γ, v, ω, λ, M₁, M₂, α , i-10 )[:]
    elseif i in [15;16;17;18;19;20;21;22;23;24]
        F[:] = Gᵢ!(F, a, w, p, θ₀, L, σ, M₁, M₂, α , i-14 )[:]
    elseif i == 25
        F[:] = E!(F,eigenpairs, M₁, M₂, α, ξ₁ᵏ ,ξ₂ᵏ, star)[:]
    elseif i in [26;27;28;29]
        F[:] =  Pᵢ!(F, p, eigenpairs, M₁, M₂, α, i-25 )[:]
    end


    return F

end

function Γᵢ!(Γ,γ,ω,M₁,M₂,α,i) 
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
    D = Derivative(1) 
    if i == 1
        Γ₁ = project(-(D*γ₁) + ω*γ₂ , space(Γ))
        Γ[:] =  Γ₁[:] 
    elseif i == 2
        Γ₂ = project((D*γ₂), space(Γ))
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        Γ₂ = Γ₂ - project(α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(Γ))
                    end
                end
            end
        end
        Γ[:] = Γ₂[:] 
    elseif i == 3
        Γ₃ = project(-(D*γ₃) + ω*γ₄, space(Γ))
        Γ[:] =  Γ₃[:] 
    elseif i == 4
        Γ₄ = project((D*γ₄), space(Γ))
        for n  in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[n,m,ℓ] != 0
                        Γ₄ = Γ₄ - project(α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(Γ))
                    end
                end
            end
        end
        Γ[:] =  Γ₄[:] 
    end

    return Γ
end

function Vᵢ!( V, ω, v, γ, M₁, M₂,α, λ, κ, lengthᵥ,i)
    γ_ext =  γ_extened_orientable!(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    v_κ = component(v,κ)

    γ₁ = component(γ_ext,1)
    γ₃ = component(γ_ext,3)

    D = Derivative(1)

    if i == 1
        V[1] =  evaluate(v_κ,0) - lengthᵥ
    elseif i == 2
        V[:] = project((D*v₁) + λ*v₁ - ω*v₂,space(V))[:]
    elseif i == 3
        V₂ = project((D*v₂) + λ*v₂,space(V))
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-2 >= 0 && m-2>= 0
                            V₂ = V₂ - ω*M₁[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V))
                        elseif n-2 >= 0
                            V₂ = V₂ - ω*M₁[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ ,space(V))
                        elseif m-2 >= 0
                            V₂ = V₂ - ω*M₁[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V))
                        end
                    end
                end
            end
        end
        V[:] = V₂[:]
    elseif i == 4
        V[:] = project((D*v₃) + λ*v₃ - ω*v₄,space(V))[:]
    elseif i == 5
        V₄ = project((D*v₄) + λ*v₄ ,space(V))
        for n  in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                     if M₂[n,m,ℓ] != 0
                        if n-2 >= 0 && m-2>= 0
                            V₄ = V₄ - ω*M₂[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃  ,space(V))
                        elseif n-2 >= 0
                            V₄ = V₄ - ω*M₂[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁  ,space(V))
                        elseif m-2 >= 0
                            V₄ = V₄ - ω*M₂[n,m,ℓ] * α^ExactReal(ℓ-1)*project(ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V))
                        end
                    end
                end
            end
        end
        V = V₄[:]
    end
    return V
end

function Wᵢ!(W, w, γ, v, ω, λ, M₁, M₂,α, k )
    γ_ext =  γ_extened_orientable!(γ)
    w₁,w₂,w₃,w₄ = eachcomponent(w)
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ_ext)
    v₁,v₂,v₃,v₄ = eachcomponent(v)

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    if eltype(α) == Interval{Float64}
        M = Sequence(Fourier(0,interval(1)) ⊗ Taylor(1), interval.([0, 1]))
    else
        M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])
    end


    if k == 1
        W₁ = project(D₁₀*w₁ + λ*M*(D₀₁*w₁) - ω*w₂, space(W))
        W₁[(:,0)] = project(Sequence( space(w₁)[1], w₁[(:,0)] ) - γ₁,space(W₁)[1])[:]
        W₁[(:,1)] = project(Sequence( space(w₁)[1], w₁[(:,1)] ) - v₁,space(W₁)[1])[:]
        W[:] =  W₁[:]
    elseif k == 2
        W₂ = project(D₁₀*w₂ + λ*M*(D₀₁*w₂) , space(W))
        for i in axes(M₁,1)
            for j in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[i,j,ℓ] != 0
                        W₂ = W₂ - α^ExactReal(ℓ-1)*ω*M₁[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(W))
                    end
                end
            end
        end
        W₂[(:,0)] = project(Sequence( space(w₂)[1], w₂[(:,0)] ) - γ₂,space(W₂)[1])[:]
        W₂[(:,1)] = project(Sequence( space(w₂)[1], w₂[(:,1)] ) - v₂,space(W₂)[1])[:]
        W[:] =  W₂[:]
    elseif k == 3
        W₃ = project(D₁₀*w₃ + λ*M*(D₀₁*w₃) - ω*w₄, space(W))
        W₃[(:,0)] = project(Sequence( space(w₃)[1], w₃[(:,0)] ) - γ₃,space(W₃)[1])[:]
        W₃[(:,1)] = project(Sequence( space(w₃)[1], w₃[(:,1)] ) - v₃,space(W₃)[1])[:]
        W[:] =   W₃[:] 
    elseif k == 4
        W₄ = project(D₁₀*w₄ + λ*M*(D₀₁*w₄) , space(W))
        for i  in axes(M₂,1)
            for j in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[i,j,ℓ] != 0
                        W₄ = W₄ - α^ExactReal(ℓ-1)*ω*M₂[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(W))
                    end
                end
            end
        end
        W₄[(:,0)] = project(Sequence( space(w₄)[1], w₄[(:,0)] ) - γ₄,space(W₄)[1])[:]
        W₄[(:,1)] = project(Sequence( space(w₄)[1], w₄[(:,1)] ) - v₄,space(W₄)[1])[:]
        W[:] =   W₄[:]
    end
    return W
end

function Gᵢ!(G, a, w, p, θ, L, σ, M₁, M₂, α , k)
    a₁, a₂, a₃, a₄ = eachcomponent(a) 
    p₁, p₂, p₃, p₄ = eachcomponent(p) 
    w₁, w₂, w₃, w₄ = eachcomponent(w) 
    σ₁, σ₂ = eachcomponent(σ)
    θ₁, θ₂ = eachcomponent(θ)

    if k == 1
        if eltype(α) == Interval{Float64}
            G[1] = p₁(σ₁[1]+interval(1im)*σ₂[1] ,σ₁[1]-interval(1im)*σ₂[1]) - a₁(1)
        else
            G[1] = p₁(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₁(1)
        end
    elseif k == 2
        if eltype(α) == Interval{Float64}
            G[1] = p₂(σ₁[1]+interval(1im)*σ₂[1] ,σ₁[1]-interval(1im)*σ₂[1]) - a₂(1)
        else
            G[1] = p₂(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₂(1)
        end
    elseif k == 3
        if eltype(α) == Interval{Float64}
            G[1] = p₃(σ₁[1]+interval(1im)*σ₂[1] ,σ₁[1]-interval(1im)*σ₂[1]) - a₃(1)
        else
            G[1] = p₃(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₃(1)
        end
    elseif k == 4
        if eltype(α) == Interval{Float64}
            G[1] = p₄(σ₁[1]+interval(1im)*σ₂[1] ,σ₁[1]-interval(1im)*σ₂[1]) - a₄(1)
        else
            G[1] = p₄(σ₁[1]+1im*σ₂[1] ,σ₁[1]-1im*σ₂[1]) - a₄(1)
        end
    elseif k == 5
        G[:] = project( ExactReal.(2*(0:order(a₁))).*a₁ , space(G) )[:]
        Ψ₁ = L/ExactReal(2)*project( a₂ ,Chebyshev(order(G)+1) )
        G[1:end] = G[1:end] + Ψ₁[2:end]  - Ψ₁[0:end-2] 
        N_w = order(w₁)    
        n = Sequence(space(w₁),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        m = Sequence(space(w₁),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        if eltype(α) == Interval{Float64}
            G[0] = sum(w₁[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^ExactReal.(n[(1:N_w[1],:)])).*interval(-0.95).^ExactReal.(m[(1:N_w[1],:)])+ w₁[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^ExactReal.(abs.(n[(-N_w[1]:-1,:)]))).*interval(-0.95).^ExactReal.(m[(-N_w[1]:-1,:)]) ) + sum( w₁[(0,:)].*interval(-0.95).^ExactReal.(m[(0,:)]))
        else
            G[0] = sum(w₁[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₁[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₁[(0,:)].*(-0.95).^m[(0,:)])
        end
        G[0] = G[0]  - a₁(-1)
    elseif k == 6
        G[:] = project( ExactReal.(2*(0:order(a₂))).*a₂ , space(G) )[:]
        Ψ₂ = Sequence(Chebyshev(order(G)+1),zeros(eltype(G), order(G)+2))
        for i in axes(M₁,1)
            for j in axes(M₁,2)
                for ℓ  in axes(M₁,3)
                    if M₁[i,j,ℓ] != 0
                        Ψ₂ = Ψ₂ + L/ExactReal(2)*M₁[i,j,ℓ] * α^ExactReal(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₂))
                    end
                end
            end
        end
        G[1:end] = G[1:end] + Ψ₂[2:end]  - Ψ₂[0:end-2] 
        N_w = order(w₂)    
        n = ExactReal.(Sequence(space(w₂),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        m = ExactReal.(Sequence(space(w₂),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        if eltype(α) == Interval{Float64}
            G[0] = sum(w₂[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)] + w₂[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₂[(0,:)].*interval(-0.95).^m[(0,:)])
        else
            G[0] = sum(w₂[(1:N_w[1],:)].*((θ₁[1] + (1im)*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₂[(-N_w[1]:-1,:)].*((θ₁[1] - (1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₂[(0,:)].*(-0.95).^m[(0,:)])
        end
        G[0] =  G[0] - a₂(-1)

    elseif k == 7
        G[:] = project( ExactReal.(2*(0:order(a₃))).*a₃ , space(G) )[:]
        Ψ₃ = L/ExactReal(2)*project( a₄ ,Chebyshev(order(G)+1) )
        G[1:end] = G[1:end] + Ψ₃[2:end]  - Ψ₃[0:end-2] 
        N_w = order(w₃)    
        n = ExactReal.( Sequence(space(w₃),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        m = ExactReal.(Sequence(space(w₃),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        if eltype(α) == Interval{Float64}
            G[0] = sum(w₃[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)] + w₃[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₃[(0,:)].*interval(-0.95).^m[(0,:)])
        else
            G[0] = sum(w₃[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₃[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₃[(0,:)].*(-0.95).^m[(0,:)])
        end
        G[0] = G[0] - a₃(-1)


    elseif k == 8
        G[:] = project( ExactReal.(2*(0:order(a₄))).*a₄ , space(G) )[:]
        Ψ₄ = Sequence(Chebyshev(order(G)+1),zeros(eltype(G),order(G)+2))
        for i  in axes(M₂,1)
            for j in axes(M₂,2)
                for ℓ  in axes(M₂,3)
                    if M₂[i,j,ℓ] != 0
                        Ψ₄ = Ψ₄ + L/ExactReal(2)*M₂[i,j,ℓ] * α^ExactReal(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₄))
                    end
                end 
            end
        end
        G[1:end] = G[1:end] + Ψ₄[2:end]  - Ψ₄[0:end-2] 
        N_w = order(w₄)    
        n = ExactReal.(Sequence(space(w₄),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        m = ExactReal.(Sequence(space(w₄),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
        if eltype(α) == Interval{Float64}
            G[0] = sum(w₄[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)] + w₄[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₄[(0,:)].*interval(-0.95).^m[(0,:)])
        else
            G[0] = sum(w₄[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)] + w₄[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)] ) + sum( w₄[(0,:)].*(-0.95).^m[(0,:)])
        end
      
        G[0] = G[0] - a₄(-1)

    elseif k == 9
        G[1] = ExactReal(0.95) - σ₁[1]^ExactReal(2) - σ₂[1]^ExactReal(2)
    elseif k == 10
        G[1] = ExactReal(1) - θ₁[1]^ExactReal(2) - θ₂[1]^ExactReal(2)
    end
    return G
end

function Pᵢ!(P, p, eigenpairs, M₁, M₂, α , k)
    p₁, p₂, p₃, p₄ = eachcomponent(p)
    λ₁ = component(eigenpairs,1)[1]
    ξ₁ = component(eigenpairs,2:5)[:]
    λ₂ = component(eigenpairs,6)[1]
    ξ₂ = component(eigenpairs,7:10)[:] 

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    if eltype(α) == Interval{Float64}
        M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), interval.([0, 1]))
        M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), interval.([0, 1]))
    else
        M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
        M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])
    end


    if k == 1
        P₁ = project(λ₁*M₁₀*(D₁₀*p₁) + λ₂*M₀₁*(D₀₁*p₁) - p₂ , space(P))
        P₁[(0,0)] = p₁[(0,0)] 
        P₁[(1,0)] = p₁[(1,0)] - ξ₁[1]
        P₁[(0,1)] = p₁[(0,1)] - ξ₂[1]
        P[:] =  P₁[:]
    elseif k ==2
        P₂ = project(λ₁*M₁₀*(D₁₀*p₂) + λ₂*M₀₁*(D₀₁*p₂) , space(P))
        for i in axes(M₁,1)
            for j in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[i,j,ℓ] != 0
                        P₂ = P₂ - M₁[i,j,ℓ] * α^ExactReal(ℓ-1)*project(p₁^(i-1)*p₃^(j-1), space(P))
                    end
                end
            end
        end
        P₂[(0,0)] = p₂[(0,0)] 
        P₂[(1,0)] = p₂[(1,0)] - ξ₁[2]
        P₂[(0,1)] = p₂[(0,1)] - ξ₂[2]
        P[:] =  P₂[:]
    elseif k == 3
        P₃ = project(λ₁*M₁₀*(D₁₀*p₃) + λ₂*M₀₁*(D₀₁*p₃) - p₄ , space(P))
        P₃[(0,0)] = p₃[(0,0)] 
        P₃[(1,0)] = p₃[(1,0)] - ξ₁[3]
        P₃[(0,1)] = p₃[(0,1)] - ξ₂[3]
        P[:] =  P₃[:]
    elseif k == 4
        P₄ = project(λ₁*M₁₀*(D₁₀*p₄) + λ₂*M₀₁*(D₀₁*p₄) , space(P))
        for i  in axes(M₂,1)
            for j in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[i,j,ℓ] != 0
                        P₄ = P₄ - M₂[i,j,ℓ] * α^ExactReal(ℓ-1)*project(p₁^(i-1)*p₃^(j-1), space(P))
                    end
                end
            end
        end
        P₄[(0,0)] = p₄[(0,0)] 
        P₄[(1,0)] = p₄[(1,0)] - ξ₁[4]
        P₄[(0,1)] = p₄[(0,1)] - ξ₂[4]
        P[:] =  P₄[:]
    
    end 
    return P
end

function V2!( V, ω, x, γ, M₁, M₂,α, κ, lengthᵥ)
    λ = component(x,1)[1]
    v = component(x,2:5)
    c₀ = component(x,6)[1]
    c = component(x,7:10)

    γ_ext =  γ_extened_orientable!(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    c₁,c₂,c₃,c₄ = eachcomponent(c)

    v_κ = component(v,κ)
    c_κ = component(c,κ)


    γ₁ = component(γ_ext,1)
    γ₃ = component(γ_ext,3)

    D = Derivative(1)
    
    V₀ =  evaluate(v_κ,0) - lengthᵥ
    V₁ = project((D*v₁) + λ*v₁ - ω*v₂,space(V)[2])
    V₂ = project((D*v₂) + λ*v₂,space(V)[3])
    V₃ = project((D*v₃) + λ*v₃ - ω*v₄,space(V)[4])
    V₄ = project((D*v₄) + λ*v₄ ,space(V)[5])

    C₀ =  evaluate(c_κ,0) 
    C₁ = project(c₀*v₁ + (D*c₁) + λ*c₁ - ω*c₂,space(V)[2])
    C₂ = project(c₀*v₂ + (D*c₂) + λ*c₂,space(V)[3])
    C₃ = project(c₀*v₃ + (D*c₃) + λ*c₃ - ω*c₄,space(V)[4])
    C₄ = project(c₀*v₄ + (D*c₄) + λ*c₄ ,space(V)[5])

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n-2 >= 0 && m-2>= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[3])
                        C₂ = C₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃ ,space(V)[8])
                    elseif n-2 >= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ ,space(V)[3])
                        C₂ = C₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((n-1)*c₁^(n-2)*γ₃^(m-1)*c₁ ,space(V)[8])
                    elseif m-2 >= 0
                        V₂ = V₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[3])                      
                        C₂ = C₂ - ω*M₁[n,m,ℓ] * α^(ℓ-1)*project((m-1)*c₁^(n-1)*γ₃^(m-2)*c₃ ,space(V)[8])
                    end
                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                 if M₂[n,m,ℓ] != 0
                    if n-2 >= 0 && m-2>= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃  ,space(V)[5])
                        C₄ = C₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃  ,space(V)[10])
                    elseif n-2 >= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁  ,space(V)[5])
                        C₄ = C₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁  ,space(V)[10])
                    elseif m-2 >= 0
                        V₄ = V₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(V)[5])
                        C₄ = C₄ - ω*M₂[n,m,ℓ] * α^(ℓ-1)*project((m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃ ,space(V)[10])
                    end
                end
            end
        end
    end

    V[:] = [V₀; V₁[:] ; V₂[:] ; V₃[:] ; V₄[:] ; C₀; C₁[:] ; C₂[:] ; C₃[:] ; C₄[:]] 
    return V
end