function PAC!( X_init, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂ ,N, max_arclength)
    X = copy(X_init)
    α₀ = real(X[1])
    x₀ = component(X,2:29)
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x₀)

    display(scatter([α₀],[L], color = "red", legend = false,markersize = 3,size=(750,750) ))

    space_X = space(X)

    f₀ = zeros(ComplexF64,space_f!(x₀))
    Dₓf₀ = LinearOperator( space(x₀), space(f₀) , zeros(ComplexF64, length(f₀), length(x₀) )  )  
    Dαf₀ = zero(f₀) 
    space_f = space(f₀)
    space_F = ParameterSpace() × space(f₀)
    
    N_fft = nextpow(2, 2N + 1)
    npts = N_fft ÷ 2 + 1

    arclength = max_arclength
    arclength_grid = [0.5 * arclength - 0.5 * cospi(2j/N_fft) * arclength for j ∈ 0:npts-1]

    X_grid = Vector{typeof(X)}(undef, npts)
    Xₖ_dot_grid = Vector{typeof(X)}(undef, npts)

    direction = zero(X)
    direction[1] = -1 # Starting direction od the paramter

    X_grid[1] = X
    Dαf = DαF!( copy(Dαf₀) , component(X_grid[1],2:29)  , M₁, M₂, real(X_grid[1][1]) )
    Dₓf = DF!( copy(Dₓf₀), component(X_grid[1],2:29)  , κ, star, M₁, M₂, real(X_grid[1][1]))
    DXf = LinearOperator( space_X , space_f ,  [ Dαf[:]  Dₓf[:,:]] )


    Q, = qr(transpose(conj(Complex.(DXf[:,:]))));
    Xₖ_dot  = Q[:,end]*conj(Q[1,end])/abs(Q[1,end])
    Xₖ_dot[2:end] = x_remove_complex!(Sequence(space_x,Xₖ_dot[2:end]))[:]


    Xₖ_dot_grid[1] = Sequence(space_X,Xₖ_dot )

    if real(direction ⋅ Xₖ_dot_grid[1]) < 0 # enforce direction
        Xₖ_dot_grid[1] .*= -1
    end
    
    for i ∈ 2:npts
        δᵢ = arclength_grid[i] - arclength_grid[i-1]
    
        Xₖ = X_grid[i-1] .+ δᵢ .* Xₖ_dot_grid[i-1]

        X = Newton_F_all!( X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , Xₖ_dot_grid[i-1], 5*10^-8 , 50 )
        α = real(component(X,1)[1])
        x = component(X,2:29)
        γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
        
        x,k = Newton_F!( x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α, 5*10^-13 , 10)
        X = Sequence( ParameterSpace() × space(x) , [α;x[:]] )
        

        display(scatter!([real(X[1])],[real(L)], color = "red", legend = false,markersize = 3,size=(750,750) ))

    
        X_grid[i] = X

        Dαf = DαF!( copy(Dαf₀) , component(X_grid[i],2:29)  , M₁, M₂, real(X_grid[i][1]) )
        Dₓf = DF!( copy(Dₓf₀), component(X_grid[i],2:29)  , κ, star, M₁, M₂, real(X_grid[i][1]))
        DXf = LinearOperator( space_X , space_f ,  [ Dαf[:]  Dₓf[:,:]] )
        
        Q, = qr(transpose(conj(Complex.(DXf[:,:]))));
        Xₖ_dot  = Q[:,end]*conj(Q[1,end])/abs(Q[1,end])
        Xₖ_dot[2:end] = x_remove_complex!(Sequence(space_x,Xₖ_dot[2:end]))[:]
        Xₖ_dot_grid[i] = Sequence(space_X,Xₖ_dot )
        #=
        Q = vec(nullspace(DXf[:,:]))
        Xₖ_dot  = Q[:]*conj(Q[1])/abs(Q[1])
        Xₖ_dot_grid[i] = Sequence(space_X,Xₖ_dot )
        =#
        if real(Xₖ_dot_grid[i-1] ⋅ Xₖ_dot_grid[i]) < 0 # enforce direction
            Xₖ_dot_grid[i] .*= -1
        end
    end
    return X_grid, Xₖ_dot_grid
    
end
