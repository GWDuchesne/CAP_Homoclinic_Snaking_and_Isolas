# ====================================================
# Newton method to find numerical candidates
# ====================================================

function Newton_γ!(γ,ω,M₁,M₂,α,ϵ,k_max) 
    γₙ = γ
    F = Γ!(Sequence( space(Derivative(1)*γ)  , zeros(length(γ),1)[:] ),γₙ,ω,M₁,M₂,α) 
    k = 0
    str = String("   at k = $(k), ||F|| = $(norm(F))")
    println(str)
    DωΓ = zero( Derivative(1)*γ)
    while norm(F) > ϵ && k <= k_max
        k = k+1
        DF,DωΓ = DΓ!(LinearOperator(space(γ),space(Derivative(1)*γ),zeros(length(γ),length(γ))),DωΓ,γₙ,ω,M₁,M₂,α)  
        γₙ = γₙ -  DF\F
        F = Γ!(Sequence( space(Derivative(1)*γ)  , zeros(length(γ),1)[:] ),γₙ,ω,M₁,M₂,α)
        str = String("   at k = $(k), ||F|| = $(norm(F))")
        println(str)
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return γₙ
end
function Newton_v!(  ω, v, γ, M₁, M₂,α, λ, κ, lengthᵥ,ϵ, k_max)
    λₙ = λ
    vₙ = v
    xₙ = Sequence(ParameterSpace() × space(v),[λₙ;vₙ[:]])
    F = V!( Sequence( ParameterSpace() × space(v),zeros(ComplexF64,length(v)+1,1)[:] ) , ω, vₙ, γ, M₁, M₂,α, λₙ, κ, lengthᵥ)
    k = 0
    display(norm(F))
    while norm(F) > ϵ && k <= k_max
        DλV = Sequence( ParameterSpace() × space(v), zeros(ComplexF64,length(v)+1,1)[:])
        DᵥV = LinearOperator(space(v),ParameterSpace() × space(v),zeros(ComplexF64,length(v)+1,length(v)))
        DωV = Sequence( ParameterSpace() × space(v), zeros(ComplexF64,length(v)+1,1)[:])
        DᵧV = LinearOperator(space(γ),ParameterSpace() × space(v),zeros(ComplexF64,length(v)+1,length(γ)))
        DλF, DᵥF, DωV, DᵧV = DV!(DλV, DᵥV, DωV, DᵧV , ω, v, γ, M₁, M₂,α, λ, κ)
        DF = LinearOperator(ParameterSpace() × space(v),ParameterSpace() × space(v),[ DλF[:] DᵥF[:,:]])
        xₙ = xₙ -  DF\F
        λₙ = component(xₙ,1)[1]
        vₙ = Sequence(space(v),component(xₙ,2:5)[:])
        xₙ = Sequence(ParameterSpace() × space(v),[λₙ;vₙ[:]])
        F = V!( Sequence( ParameterSpace() × space(v),zeros(ComplexF64,length(v)+1,1)[:] ) , ω, vₙ, γ, M₁, M₂,α, λₙ, κ, lengthᵥ)
        display(norm(F))
        k = k+1
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return λₙ,vₙ
end
function Newton_E!(eigenpairs, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star) 
    eigenpairsₙ = eigenpairs
    F = E!(zero(eigenpairs),eigenpairsₙ, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star)
    k = 0
    display(norm(F))
    while norm(F) > ϵ && k <= k_max
        DF = LinearOperator( space(E),space(E) ,zeros(ComplexF64,10,10))
        DF = DE!(DF, eigenpairsₙ, M₁, M₂,α, star)
        eigenpairsₙ = eigenpairsₙ -  DF\F
        F = E!(zero(eigenpairs),eigenpairsₙ, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star)
        display(norm(F))
        k = k+1
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return eigenpairsₙ
end
function Newton_W!( w, γ_ext, v, ω, λ, M₁, M₂,α,ϵ, k_max)
    wₙ = w
    F = W!(zero(w), wₙ, γ_ext, v, ω, λ, M₁, M₂,α)
    k = 0
    display(norm(F))
    while norm(F) > ϵ && k <= k_max
        DωW = zero(w)
        DλW = zero(w)
        DᵥW = LinearOperator( space(v) , space(w) , zeros(ComplexF64,length(w),length(v)))
        DwW = LinearOperator( space(w) , space(w) , zeros(ComplexF64,length(w),length(w)))
        DᵧW = LinearOperator( space(γ) , space(w) , zeros(ComplexF64,length(w),length(γ)))
        DωW, DλW, DᵥW, DF, DᵧW  = DW!(DωW, DλW, DᵥW, DwW, DᵧW, wₙ, γ, v, ω, λ, M₁, M₂,α)
        wₙ = wₙ -  DF\F
        F = W!(zero(w), wₙ, γ_ext, v, ω, λ, M₁, M₂,α)
        display(norm(F))
        k = k+1
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return wₙ
end
function Newton_P!( p, eigenpairs, M₁, M₂, α)
    pₙ = p
    F = P!(zero(p), pₙ, eigenpairs, M₁, M₂, α)
    k = 0
    display(norm(F))
    while norm(F) > ϵ && k <= k_max
        DₚP = LinearOperator( space(p), space(p), zeros(ComplexF64,length(p),length(p)) )
        DeigenpairsP = LinearOperator( space(eigenpairs), space(p), zeros(ComplexF64,length(p),length(eigenpairs)) )
        
        DeigenpairsP ,DF = DP!(DeigenpairsP, DₚP,pₙ, eigenpairs , M₁, M₂, α)
        pₙ = pₙ -  DF\F
        F = P!(zero(p), pₙ, eigenpairs, M₁, M₂, α)
        display(norm(F))
        k = k+1
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return pₙ
end

function Newton_G!(a, w, p, θ₀, L, σ, M₁, M₂, α,ϵ,k_max)
    aₙ = a
    G₀ = zeros( ComplexF64, ParameterSpace()  × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace()× ParameterSpace()) 
    G = G!(copy(G₀), aₙ, w, p, θ₀, L, σ, M₁, M₂, α)
    F = component(G,5:8)
    k = 0
    display(norm(F))
    while norm(F) > ϵ && k <= k_max
        DₐG = LinearOperator( space(a),space(G),zeros(length(G),length(a)) )  
        Dθ₀G = LinearOperator(  space(σ),space(G₀), zeros(Complex,length(G₀),2) )    
        DLG = zero(G₀)
        DσG = LinearOperator(  space(σ),space(G₀), zeros(Complex,length(G₀),2) )
        DₚG = LinearOperator( space(p),space(G₀), zeros(Complex,length(G₀),length(p)) )
        DwG = LinearOperator( space(w),space(G₀), zeros(Complex,length(G₀),length(w)) )
        DₐG, Dθ₀G, DLG, DσG, DₚG, DwG =  DG!(DₐG, Dθ₀G, DLG, DσG, DₚG, DwG, aₙ, w, p, θ₀, L, σ, M₁, M₂, α)
        DF = component(DₐG,5:8,:)
        aₙ = aₙ -  real(DF\F)
        G = G!(copy(G₀), aₙ, w, p, θ₀, L, σ, M₁, M₂, α)
        F = component(G,5:8)
        display(norm(F))
        k = k+1
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return aₙ
end

function Newton_F!( x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α, ϵ , k_max)
    xₙ = x
    F₀ = zeros(ComplexF64,space_f!(x))
    F = F!( F₀ , xₙ, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
    k = 0
    str = String("   at k = $(k), ||F|| = $(norm(F))")
    println(str)
    while norm(F) > ϵ && k <= k_max
        DF = LinearOperator( space(x), space(F₀) , zeros(ComplexF64, length(F₀), length(x) )  )        
        DF = DF!(DF, xₙ, κ, star, M₁, M₂, α)
        xₙ = xₙ -  x_remove_complex!(DF\F)
        F = F!( F₀ , xₙ, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α)
        k = k+1
        str = String("   at k = $(k), ||F|| = $(norm(F))")
        println(str)
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return xₙ,k
end

function Newton_F_all!( X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot, ϵ, k_max )
    
    Xₙ = X
    αₙ = real(X[1])
    xₙ = x_remove_complex!(component(X,2:29))
    F₀ = zeros(ComplexF64,ParameterSpace() × space_f!(xₙ))
    F = F_all!(copy(F₀), Xₙ, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot )
    k = 0
    norm_F = 0
    println("   At k = $(k), we have ||F|| = $(norm(F))")
    while norm(F) > ϵ && k <= k_max
        DF = LinearOperator( space(X) , space(F₀) , zeros(ComplexF64,length(X),length(F)))
        DF = DF_all!(DF, X, κ, star, lengthᵥ, M₁, M₂ , Xₖ_dot )
        Xₙ = Xₙ - DF\F

        αₙ = real(component(Xₙ,1)[1])
        xₙ = x_remove_complex!(component(Xₙ,2:29))
        Xₙ[:] = Sequence( space(Xₙ), [αₙ;xₙ[:]])[:]

        F = F_all!(copy(F₀), Xₙ, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot )
        k = k+1
        println("   At k = $(k), we have ||F|| = $(norm(F))")
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return Xₙ
end

function Newton_EP!( p, eigenpairs, M₁, M₂, α, ξ₁ᵏ ,ξ₂ᵏ, star , ϵ, k_max)
    space_x = ParameterSpace()^10 × space(p)

    eigenpairsₙ = eigenpairs
    pₙ = p
    xₙ = Sequence(space_x,[eigenpairsₙ[:]; pₙ[:]])

    E = E!(zero(eigenpairs),eigenpairsₙ, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star)
    P = P!(zero(p), pₙ, eigenpairsₙ, M₁, M₂, α)
    F = Sequence(space_x,[E[:];P[:]])

    k = 0
    str = String("   at k = $(k), ||F|| = $(norm(F))")
    println(str)
    while norm(F) > ϵ && k <= k_max
        DF = LinearOperator( space_x, space_x, zeros(ComplexF64,length(F),length(F)) )
        DF = DEP(  DF, α, eigenpairsₙ, pₙ , star, M₁, M₂) 
        xₙ = xₙ -  DF\F
        eigenpairsₙ = component(xₙ,1)
        pₙ = component(xₙ, 2:5)

        E = E!(zero(eigenpairs),eigenpairsₙ, M₁, M₂,α, ξ₁ᵏ ,ξ₂ᵏ, star)
        P = P!(zero(p), pₙ, eigenpairsₙ, M₁, M₂, α)
        F = Sequence(space_x,[E[:];P[:]])
        k = k+1
        str = String("   at k = $(k), ||F|| = $(norm(F))")
        println(str)
        if k > k_max
            error("Reached the maximum numbers of steps")
        end
    end

    return eigenpairsₙ, pₙ
end