#---------------------------------------------------------------------
# Solutions Finder Helper
#---------------------------------------------------------------------
    # The goal of this file is to help users to find a numerical 
    # approximation of a heteroclinic orbit connecting a periodic solution
    # to the equilibrium at the origin

#---------------------------------------------------------------------
# Packages 
#---------------------------------------------------------------------
    # Add the Following Packages if not already installed:
    #  - RadiiPolynomial
    #  - LinearAlgebra
    #  - JLD2
    #  - Plots
    #  - DifferentialEquations
    #  - Colors
    #  - ColorSchemes
    #  - Polynomials

    using RadiiPolynomial, LinearAlgebra, JLD2, Plots, Colors, ColorSchemes,DifferentialEquations, Polynomials

#---------------------------------------------------------------------
# Functions Files
#---------------------------------------------------------------------
    # This add functions from other files
    include("Other_functions.jl")
    include("Map_functions_Orientable.jl")
    include("Map_derivative_Orientable.jl")
    include("Newton_functions.jl")
    include("plots and visuals.jl")

#---------------------------------------------------------------------
# 0. Setup of System
#---------------------------------------------------------------------
    printstyled("0. Setup of System\n"; color = :blue)
    # Defining Matrix M₁ and M₂ such that
    #   F₁(u,v) = ∑ᵢⱼₖ M₁[i,j,k]  αᵏ uⁱ vʲ 
    #   F₂(u,v) = ∑ᵢⱼₖ M₂[i,j,k]  αᵏ uⁱ vʲ 

    # Example for Grey-Scott equation

    m = 0.5
    d = 9
    a = 1.015
    α = ((a/m) + sqrt( (a/m )^2 -4 ) )/2
    

    M₁ = zeros(Float64,2,3,4)
    M₁[:,:,1] = [   0.0 0.0 m;
                    0.0 0.0 0.0;]

    M₁[:,:,2] = [   0.0 2*m 0.0;
                    1.0 0.0 1.0;]

    M₁[:,:,3] = [   0.0 0.0 0.0;
                    0.0 2.0 0.0;]

    M₁[:,:,4] = [   0.0 0.0 0.0;
                    1.0 0.0 0.0;]


    M₁ = M₁ ./ d

    M₂ = zeros(Float64,2,3,4)
    M₂[:,:,1] = [   0.0 0.0 -m;
                    0.0 0.0 0.0;]

    M₂[:,:,2] = [   0.0 -m  0.0;
                    0.0 0.0 -1.;]

    M₂[:,:,3] = [   0.0 0.0 0.0;
                    0.0 -2. 0.0;]

    M₂[:,:,4] = [   0.0 0.0 0.0;
                    -1. 0.0 0.0;]



    # Parameter α
    strα = String("   α = $(α)")
    strf₁= String("   uₓₓ =")
    strf₂= String("   vₓₓ =")
    println(strα)
    println(strf₁*str_system( M₁ ))
    println(strf₂*str_system( M₂ ))


#---------------------------------------------------------------------
# 1. Local Stable Manifold of Equilibrium at the Origin
#---------------------------------------------------------------------
    printstyled("\n1. Local Stable Manifold of Equilibrium at the Origin\n"; color = :blue)
    # In this section we develop the tool to find a numerical 
    # approximation of the Local Stable MAnifold of the equilibrium at 0

    # Defining Spaces of Coefficients pᵢ and Eigenvalues/Eigenvectors
    Nₚ = [15 ; 15] # Change the number of coefficients as needed to reach ε machine
    space_pᵢ = Taylor( Nₚ[1] ) ⊗ Taylor( Nₚ[2] ) 
    space_p = space_pᵢ × space_pᵢ × space_pᵢ × space_pᵢ 
    space_eigenpairs = ParameterSpace()^10

    # Computation of DΨ(0) and eigenpair
     # Increase/Decrease this value for biggger/smaller Local manifold
    function DΨ₀!(DΨ₀, M₁, M₂ , α )
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
        DΨ₀ = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
        for  ℓ in axes(M₁,3)
            DΨ₀ = DΨ₀ + [0 0 0 0; α^(ℓ-1)*M₁_temp[2,1,ℓ] 0 α^(ℓ-1)*M₁_temp[1,2,ℓ] 0; 0 0 0 0; 0 0 0 0]
        end
        for  ℓ in axes(M₂,3)
            DΨ₀ = DΨ₀ + [0 0 0 0; 0 0 0 0; 0 0 0 0; α^(ℓ-1)*M₂_temp[2,1,ℓ] 0 α^(ℓ-1)*M₂_temp[1,2,ℓ] 0]
        end
        return DΨ₀
    end
    DΨ₀ = zeros( 4 , 4 )
    DΨ₀ = DΨ₀!(DΨ₀, M₁, M₂ , α )
    eigva = eigvals(DΨ₀)
    eigve = eigvecs(DΨ₀) 
    λ₁ =  eigva[findall( <(0) , real.(eigva)  )[1]]
    λ₂ =  eigva[findall( <(0) , real.(eigva)  )[2]]
    ξ₁ = eigve[:,findall( <(0) , real.(eigva)  )[1]] ./ 10
    ξ₂ = eigve[:,findall( <(0) , real.(eigva)  )[2]] ./ 10
    star = [findmax( abs.(ξ₁))[2] ; findmax( abs.(ξ₂))[2]] # Fixing biggest component uniqueness of eigenpair
    ξ₁ᵏ = ξ₁[star[1]]
    ξ₂ᵏ = ξ₂[star[2]]
    p = zeros(ComplexF64,space_p)
    eigenpairs = Sequence(space_eigenpairs, [  λ₁;ξ₁;λ₂;ξ₂ ] )
    println("   Newton's Method to refine Numerical Approximation:")
    eigenpairs, p = Newton_EP!( p, eigenpairs, M₁, M₂, α , ξ₁ᵏ ,ξ₂ᵏ, star, 10^-14, 25)
#---------------------------------------------------------------------
# 2. Periodic Solution
#---------------------------------------------------------------------
    printstyled("\n2. Periodic Solution\n"; color = :blue)
    # In this section we develop the tool to find a numerical 
    # approximation of the Periodic Solution

    # Shooting 
    function GS!(du, u, p, t)
        du[1] = u[2]
        du[2] = 1/9 * α * u[3]  + 1/18* u[3]^2  + 1/9 * α*u[1]  + 1/9* α^3 * u[1]  + 2/9*α^2 * u[1]*u[3]  + 1/9 * α * u[1] * u[3]^2
        du[3] = u[4]
        du[4] = -0.5 *α *u[3]  -0.5 *u[3]^2  - α^3 * u[1]  -2.0 *α^2 * u[1] * u[3]  - α *u[1] *u[3]^2
        nothing
    end

    function shooting(PP, p_graph, theta_p, tspan, ℓ)
        for i = 1:Nₜ
            σ₁_grpah = 0.95*(cos( theta_p[i] ) + 1im*sin( theta_p[i] ))
            σ₂_grpah = 0.95*(cos( theta_p[i] ) - 1im*sin( theta_p[i] ))
            p_graph[i,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
            u₀ = real( p(σ₁_grpah,σ₂_grpah) )[:]
            prob = ODEProblem(GS!, u₀, tspan)
            sol = solve(prob, Tsit5(),abstol=1e-14, reltol=1e-14)
            u₁ = zeros(Float64, length(sol.u) )
            u₂ = zeros(Float64, length(sol.u) )
            u₃ = zeros(Float64, length(sol.u) )
            u₄ = zeros(Float64, length(sol.u) )
            for j = 1:length(sol.u)
                u₁[j] = (sol.u)[j][1]
                u₂[j] = (sol.u)[j][2]
                u₃[j] = (sol.u)[j][3]
                u₄[j] = (sol.u)[j][4]
            end
            if maximum( abs.(u₁ ) ) > 0.3 || maximum( abs.(u₃  ) ) > 1
                println("When k = $(i): Bad at θ = $(theta_p[i])")
            else
                println("When k = $(i): Good")
            end
            if ℓ == 1
                u_plot = u₁ 
            elseif ℓ == 2
                u_plot = u₂
            elseif ℓ == 3
                u_plot = u₃
            else
                u_plot = u₄
            end
            #plot3d!( u₁ , u₂, u₃, legend = false) 
            plot!( sol.t , u_plot , legend = false)
            #plot!( u₁ , u₃, legend = false)
        end
        plot(PP)
    end

    # Shooting Backward From Stable Manifold 
    Nₜ = 2
    tspan = (0, -150) 
    # theta_p = range( 3.668045997617721 ,3.6719106112684923,Nₜ) # -38 going nowhere ?
    # theta_p = range(1.7053308815492894,2.6104524054262943,Nₜ) # -24 Alternate path
    # theta_p = range(3.1929620225557054, 3.3503380847018542,Nₜ)# -38
    # theta_p = range(3.4229338135790783, 3.475175215913335,Nₜ) # -54 redo for other branch
    #theta_p = range(3.445387337391182, 3.447119196783625,Nₜ) # -70 redo for other branch
    theta_p = range( 3.448202571464514, 3.448202571465317,Nₜ) # -70 redo for other branch
     
    p_graph = zeros(Nₜ,4)
    PP = plot(size = (750,750)) #scatter( [0] , [0], [0], legend = false) 

    shooting(PP, p_graph, theta_p, tspan,1)
   
    σ = [ sqrt(0.95)*cos( sum(theta_p)/Nₜ ) ; sqrt(0.95)*sin( sum(theta_p)/Nₜ )]
    p_at_σ = real( p(σ[1] +1im*σ[2],σ[1] -1im*σ[2]) )[:]
    u₀ = p_at_σ 
    prob = ODEProblem(GS!, u₀, tspan)
    sol = solve(prob, Tsit5(),abstol=1e-14, reltol=1e-14)
    u₁ = zeros(Float64, length(sol.u) )
    u₂ = zeros(Float64, length(sol.u) )
    u₃ = zeros(Float64, length(sol.u) )
    u₄ = zeros(Float64, length(sol.u) )
    for j = 1:length(sol.u)
        u₁[j] = (sol.u)[j][1]
        u₂[j] = (sol.u)[j][2]
        u₃[j] = (sol.u)[j][3]
        u₄[j] = (sol.u)[j][4]
    end
    plot(sol.t[end-1500:end] ,u₁[end-1500:end] ,legend = false)
    ~,k₁ =  findmin(u₁[end-1200:end])
    t₁ = (sol.t[end-1200:end])[k₁]
    ~,k₂ =  findmin(u₁[end-1500:end-1000])
    t₂ = (sol.t[end-1500:end-1000])[k₂]

    # Period, spaces and coefficients
    T = abs.( t₁ - t₂ ) # Approximation of the Period as starting point.
                        # For a Orientable solution use 2T
    ω = T / 2π # Rescaling
    Nᵧ = 20
    space_γᵢ_cos = CosFourier( Nᵧ , 1)
    space_γᵢ_sin = SinFourier( Nᵧ , 1)
    space_γ = space_γᵢ_cos × space_γᵢ_sin × space_γᵢ_cos × space_γᵢ_sin  
    space_Dγ =  space_γᵢ_sin × space_γᵢ_cos × space_γᵢ_sin  × space_γᵢ_cos

    # Data Swift–Hohenberg equation
    γ = zeros(space_γ)
    component(γ,1)[0:3] = [0.19668646619787655;
                            -0.1417409989122467;
                            0.015759794739307234;
                            0.00037594703167210473] 
    component(γ,2)[1:3] = [  0.203040263145214;
                            -0.04515098518480985;
                            -0.001615602785894526]  
    component(γ,3)[0:3] = [ -0.038167090252435636;
                            0.024335774508661798;
                            -0.0020978521467080423;
                            -3.928702291951463e-5] 
    component(γ,4)[1:3] = [ -0.03486035866828054;
                            0.006010236349061742;
                            0.00016883288955883935]  

    # Data Swift–Hohenberg equation
    γ = Newton_γ!(γ,ω,M₁,M₂,α,10^-14,25) 
    if norm(γ) < 10^-12
        error("Periodic Solution converged to the trivial solution. Try to use different initial conditions for Newton's Method")
    end

    plot(range(0,10*2*pi,1000), component(γ,1).(range(0,10*2*pi,1000)) .+ m/α , ylimit = [0,1.75] , legend = false, size = (1500,250))
    plot!(range(0,10*2*pi,1000), component(γ,3).(range(0,10*2*pi,1000)) .+ α )

#---------------------------------------------------------------------
# 3. Bundle and Local Unstable Manifold of Periodic Solution
#---------------------------------------------------------------------
    printstyled("\n3. Bundle and Local Unstable Manifold of Periodic Solution\n"; color = :blue)
    # In this section we develop the tool to find a numerical 
    # approximation of the Bundle and Local Unstable Manifold 
    # of Periodic Solution.

    lengthᵥ = 0.01 # Fixing biggest component uniqueness of eigenpair
    κ = 1 # Fixing component index

    # Spaces and coefficients
    Nᵥ = Nᵧ
    space_vᵢ = Fourier(Nᵥ,1)
    space_v = space_vᵢ × space_vᵢ × space_vᵢ × space_vᵢ 

    N_w = [Nᵥ , 12]
    space_wᵢ = Fourier(N_w[1],1) ⊗ Taylor(N_w[2])
    space_w = space_wᵢ × space_wᵢ × space_wᵢ × space_wᵢ 
    v = zeros(ComplexF64,space_v)
    λ = 0.6573892123582663
    component(v,1)[0] = 1
    component(v,2)[0] = 1
    component(v,3)[0] = 1
    component(v,4)[0] = 1

    w = zeros(ComplexF64, space_w)
    γ_ext =  γ_extened_orientable!(γ)
    component(w,1)[(:,0)] = component(γ_ext,1)[:]
    component(w,2)[(:,0)] = component(γ_ext,2)[:]
    component(w,3)[(:,0)] = component(γ_ext,3)[:]
    component(w,4)[(:,0)] = component(γ_ext,4)[:]

    component(w,1)[(:,1)] = component(v,1)[:]
    component(w,2)[(:,1)] = component(v,2)[:]
    component(w,3)[(:,1)] = component(v,3)[:]
    component(w,4)[(:,1)] = component(v,4)[:]

    # Functions for Newton's Method
    function ΓVW!( X, M₁, M₂, α,  κ, lengthᵥ ,space_out )
        γ =  component(X,1:4)
        λ =  real(component(X,5)[1])
        v =  component(X,6:9)
        w =  component(X,10:13)
        F = Sequence( space_out ,  [   Γ!(zero(Derivative(1)*γ ),γ,ω,M₁,M₂,α)[:] ;  V!( zeros(eltype(v),ParameterSpace()× space(v)), ω, v, γ, M₁, M₂, α, λ, κ, lengthᵥ)[:] ; W!( zeros(ComplexF64,space(w)), w, γ, v, ω, λ, M₁, M₂,α)[:] ]  )
        return F
    end
    function DΓVW!( DF, X, M₁, M₂, α,  κ )
        γ =  component(X,1:4)
        λ =  real(component(X,5)[1])
        v =  component(X,6:9)
        w =  component(X,10:13)

        D = Derivative(1)    
        DωΓ = zeros(  eltype(v),space(D*γ))
        DᵧΓ = zeros(eltype(v), space(γ), space(D*γ))
        DλΓ = zeros( eltype(v), space(D*γ))
        DᵥΓ = zeros( eltype(v), space(v), space(D*γ))
        DwΓ = zeros( eltype(v), space(w), space(D*γ))

        DωV = zeros(eltype(v), ParameterSpace() ×  space(v))
        DᵧV = zeros(eltype(v), space(γ), ParameterSpace() ×  space(v))
        DλV = zeros(eltype(v), ParameterSpace() × space(v))
        DᵥV = zeros(eltype(v), space(v), ParameterSpace() × space(v))
        DwV = zeros( eltype(v), space(w), ParameterSpace() × space(v))

        DωW = zeros(eltype(v), space(w))
        DᵧW = zeros(eltype(v), space(γ), space(w))
        DλW = zeros(eltype(v), space(w))
        DᵥW = zeros(eltype(v), space(v), space(w))
        DwW = zeros( eltype(v), space(w), space(w))

        DλV, DᵥV, DωV, DᵧV = DV!( DλV, DᵥV,DωV,DᵧV , ω, v, γ, M₁, M₂,α, λ, κ)
        DᵧΓ,DωΓ = DΓ!(DᵧΓ,DωΓ,γ,ω,M₁,M₂,α) 
        DωW, DλW, DᵥW, DwW, DᵧW = DW!(DωW, DλW, DᵥW, DwW, DᵧW , w, γ, v, ω, λ, M₁, M₂,α)

        DF[:,:] = [ DᵧΓ[:,:] DλΓ[:] DᵥΓ[:,:] DwΓ[:,:]; 
                    DᵧV[:,:] DλV[:] DᵥV[:,:] DwV[:,:];
                    DᵧW[:,:] DλW[:] DᵥW[:,:] DwW[:,:]]
        return DF
    end
    function Newton_ΓVW!( X, M₁, M₂, α,  κ, lengthᵥ ,space_out  )
        Xₙ = X
        F = ΓVW!( Xₙ, M₁, M₂, α,  κ, lengthᵥ ,space_out )
        k = 0
        str = String("   at k = $(k), ||F|| = $(norm(F))")
        println(str)
        while norm(F) > 10^-10 && k <= 25
            DF =  DΓVW!( zeros(ComplexF64, space(X), space_out), Xₙ, M₁, M₂, α,  κ )
            Xₙ = Xₙ - DF\F    
            F = ΓVW!( Xₙ, M₁, M₂, α,  κ, lengthᵥ ,space_out )
            k = k+1
            str = String("   at k = $(k), ||F|| = $(norm(F))")
            println(str)
            if k > 25
                error("Reached the maximum numbers of steps")
            end
        end
        return Xₙ
    end

    space_X =  space_γ × ParameterSpace() × space_v × space_w 
    space_out = space_Dγ × ParameterSpace() × space_v × space_w 
    X = Sequence(space_X, [ γ[:] ; λ ; v[:]; w[:] ] )
    println("   Newton's Method to refine Numerical Approximation:")
    X = Newton_ΓVW!( X, M₁, M₂, α,  κ, lengthᵥ ,space_out  )
    γ =  real.(component(X,1:4))
    λ =  real(component(X,5)[1])
    v =  component(X,6:9)
    w =  component(X,10:13)


#---------------------------------------------------------------------
# 4. Connecting Orbit
#---------------------------------------------------------------------
    printstyled("\n4. Connecting Orbit\n"; color = :blue)
    # In this section we develop the tool to find a numerical 
    # approximation of the Connecting Orbit between the two
    # manifolds found before. To do so, we will shoot 
    # to approximate the Chebyshev coefficients of the orbit on
    # a range of L and θ. 

    function shooting_from_periodic(PP, p_graph, theta_p, tspan, ℓ)
        for i = 1:Nₜ
            p_graph[i,:] = real( w(theta_p[i],-0.95) )[:]
            u₀ = real( w(theta_p[i],-0.95) )[:]
            prob = ODEProblem(GS!, u₀, tspan)
            sol = solve(prob, Tsit5(),abstol=1e-14, reltol=1e-14)
            u₁ = zeros(Float64, length(sol.u) )
            u₂ = zeros(Float64, length(sol.u) )
            u₃ = zeros(Float64, length(sol.u) )
            u₄ = zeros(Float64, length(sol.u) )
            for j = 1:length(sol.u)
                u₁[j] = (sol.u)[j][1]
                u₂[j] = (sol.u)[j][2]
                u₃[j] = (sol.u)[j][3]
                u₄[j] = (sol.u)[j][4]
            end
            if  abs.(u₁[end] )  > 0.1 || maximum( abs.(u₃  ) ) > 1
                println("When k = $(i): Bad at θ = $(theta_p[i])")
            else
                println("When k = $(i): Good")
            end
            if ℓ == 1
                u_plot = u₁ 
            elseif ℓ == 2
                u_plot = u₂
            elseif ℓ == 3
                u_plot = u₃
            else
                u_plot = u₄
            end
            #plot3d!( u₁ , u₂, u₃, legend = false) 
            plot!( sol.t , u_plot , legend = false)
            #plot!( u₁ , u₃, legend = false)
        end
        plot(PP)
    end

    # Shooting Backward From Stable Manifold 
    Nₜ = 100
    tspan = (0, 22) 
    theta_p = range( 6.240233219068064, 7.235183080994675, Nₜ)   

    p_graph = zeros(Nₜ,4)
    PP = plot(size = (1000,750)) #scatter( [0] , [0], [0], legend = false) 
    shooting_from_periodic(PP, p_graph, theta_p, tspan,1)
   
   
    plot!(range(  -2*pi , pi ,100 ) .-pi, component( γ , 1).(range(  -2*pi , pi ,100 )) )



    L = 10 # Approximation of length of L
    N_ang = 4
    angy = exp.( 1im*pi .* range(0.8532720787527833,1.0084124567078347,10) )

    θ = Sequence( ParameterSpace() × ParameterSpace(),[real(angy[N_ang]) ; imag(angy[N_ang])] )    # Approximation of length of angle of exit
    σ = Sequence( ParameterSpace() × ParameterSpace(),[ -0.66013241921651 ; -0.7170949651889614] ) # Approximation of length of angle of entering
    # Better approximation can be obtain by numerical shooting

    # Spaces and Variables
    Nₐ = 100
    space_aᵢ = Chebyshev(Nₐ)
    space_a = space_aᵢ × space_aᵢ × space_aᵢ × space_aᵢ 

    function G_a!( a, w, p, θ, L, σ, M₁, M₂, α)
        space_X = ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace() × ParameterSpace() 
        return component( G!(zeros(ComplexF64,space_X), a, w, p, θ, L, σ, M₁, M₂, α) , 5:8)
    end

    function newtonG_a!( a, w, p, θ, L, σ, M₁, M₂, α)
        aₙ = a
        space_X = ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space(a) × ParameterSpace() × ParameterSpace()
        F = G_a!( aₙ, w, p, θ, L, σ, M₁, M₂, α)
        k = 0
        while norm(F) > 10^-10 && k <= 25
            DₐG = zeros( ComplexF64,space(a),space_X)
            DθG = zeros( ComplexF64,ParameterSpace() × ParameterSpace(), space_X)
            DLG = zeros( ComplexF64, space_X)
            DσG = zeros( ComplexF64, ParameterSpace() × ParameterSpace(),  space_X)
            DₚG = zeros( ComplexF64,space(p),space_X)
            DwG = zeros( ComplexF64,space(w),space_X)

            DₐG, DθG, DLG, DσG, DₚG, DwG = DG!(DₐG, DθG, DLG, DσG, DₚG,DwG, aₙ, w, p, θ, L, σ, M₁, M₂, α)

            DF = component(DₐG,5:8,:)
            aₙ = aₙ - DF\F    
            F = G_a!( aₙ, w, p, θ, L, σ, M₁, M₂, α)
            k = k+1
            if k > 25
                error("Reached the maximum numbers of steps")
            end
        end
        return aₙ
    end

    space_x = space_γ × ParameterSpace() × space_v × space_w × ParameterSpace() × ParameterSpace()× ParameterSpace() × ParameterSpace()× ParameterSpace() × space_a× ParameterSpace()× ParameterSpace()^10 × space_p
    x = zeros(ComplexF64, space_x)
    a = zeros(space_a)
    component(a,1)[0] = 0
    component(a,2)[0] = 0
    component(a,3)[0] = 0
    component(a,4)[0] = 0

    # Approximation of the connecting orbit
    a = real.(newtonG_a!( a, w, p, θ, L, σ, M₁, M₂, α)) 

    # Newton's Method 
    x[:] = [γ[:];λ;v[:];w[:];L;θ[:];σ[:];a[:];ω;eigenpairs[:];p[:]]
    println("   Newton's Method to refine Numerical Approximation:")
    x,k = Newton_F!( x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α, 10^-14 , 25)

    
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)

#---------------------------------------------------------------------
# 5. Plotting
#---------------------------------------------------------------------
    printstyled("\n5. Plotting\n"; color = :blue)
    println("   Plot of the Heteroclinic Orbit ")

    Nₜ = 100
    t = range( -1,1,Nₜ)
    r = range( 0,1,Nₜ) 
    theta_p = range( 0,2*pi,Nₜ)
    p_graph = zeros(Nₜ,Nₜ,4)
    γ_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    w_graph = zeros(Nₜ,Nₜ,4)

    for i = 1:Nₜ
        γ_graph[i,:] = γ(theta_p[i] )[:]
        a_graph[i,:] = a(t[i] )[:]
        for j = 1:Nₜ
            σ₁_grpah = r[j]*(cos( theta_p[i] ) + 1im*sin( theta_p[i] ))
            σ₂_grpah =  r[j]*(cos( theta_p[i] ) - 1im*sin( theta_p[i] ))
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
            w_graph[i,j,:] = real( w(theta_p[i],t[j]) )[:]
        end
    end
    
    P = plot3d( w_graph[:,:,1][:],w_graph[:,:,2][:],w_graph[:,:,3][:],linecolor = "red",alpha = 0.1 )
    plot3d!( p_graph[:,:,1][:],p_graph[:,:,2][:],p_graph[:,:,3][:],linecolor = "cyan",alpha = 0.25 , legend = false ,size=(750,750), camera=(60,30))
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")
    plot3d!( γ_graph[:,1][:],γ_graph[:,2][:],γ_graph[:,3][:],linecolor = "red")
    plot3d!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:],linecolor = "purple")
    display(plot(P))
    
#---------------------------------------------------------------------
# 6. Chebyshev Interpolation of first Segment
#---------------------------------------------------------------------
    printstyled("\n6. Chebyshev Interpolation of first Segments\n"; color = :blue)
    # In this section we compute the first segments of the Chebychev 
    # expension of the parameter with an pseudo-arclength continuation

    N = 15 # Nymber of Chebyshev coefficients
    max_arclength = 0.25 # length of the segment
    direct = -1 # Starting direction od the parameter

    # Space and Variable
    space_X = ParameterSpace() × space_x
    X = Sequence(space_X, [α;x[:]])
    str_color = "red"
    #Pseudo-Arclength Continuation (First segment)
    println("   First segment ")
    function PAC!( X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂ ,N, max_arclength,direct,str_color)
        space_X = space(X)
        space_F =  space_F!(X)
        X₀ = copy(X)
        N_fft = nextpow(2, 2N + 1)
        npts = N_fft÷2 + 1
        α = real(X[1])
        x = component(X,2:29)

        α₀ = copy(α)
        x₀ = copy(x)
        space_x = space(x)
        γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x₀)
        P2 = scatter([α₀],[L], color = str_color, legend = false,markersize = 3,size=(750,750) )
        println("   Next point at α = $(α₀) ")
        space_f = space_f!(x₀)
        f₀ = zeros(ComplexF64, space_f)
        Dₓf₀ = zeros( ComplexF64 , space_x , space_f  )  
        Dαf₀ = zeros(ComplexF64, ParameterSpace(), space_f)

        arclength = max_arclength
        arclength_grid = [0.5 * arclength - 0.5 * cospi(2j/N_fft) * arclength for j ∈ 0:npts-1]

        X_grid = Vector{typeof(X)}(undef, npts)
        η_grid = Vector{typeof(X)}(undef, npts)

        direction = zero(X)
        direction[1] = direct # Starting direction od the parameter

        X_grid[1] = X
        Dαf = DαF!( copy(Dαf₀) , component(X_grid[1],2:29)  , M₁, M₂, real(X_grid[1][1]) )
        Dₓf = DF!( copy(Dₓf₀), component(X_grid[1],2:29)  , κ, star, M₁, M₂, real(X_grid[1][1]))
        DXf = LinearOperator( space_X , space_f ,  [ Dαf[:,:]  Dₓf[:,:]] )

        Q, = qr(transpose(conj(Complex.(DXf[:,:]))));
        η  = Q[:,end]*conj(Q[1,end])/abs(Q[1,end])
        η[2:end] = x_remove_complex!(Sequence(space_x,η[2:end]))[:]

        η_grid[1] = Sequence(space_X, η )

        if real(direction ⋅ η_grid[1]) < 0 # enforce direction
            η_grid[1] .*= -1
        end
        for i ∈ 2:npts
            δᵢ = arclength_grid[i] - arclength_grid[i-1]
        
            Xₖ = X_grid[i-1] .+ δᵢ .* η_grid[i-1]

            #test = Newton_F_all!( copy(X₀), κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , η_grid[i-1], 5*10^-8 , 50 )
            X₀ = Newton_F_all!( X₀, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Xₖ , η_grid[i-1], 5*10^-15 , 100 )
            α = real(component(X₀,1)[1])
            x = component(X₀,2:29)
            γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
            
            #x,k = Newton_F!( x, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, α, 5*10^-13 , 10)
            #X₀ = Sequence( ParameterSpace() × space(x) , [α;x[:]] )
            

            scatter!([real(X₀[1])],[L], color = str_color, legend = false,markersize = 3,size=(750,750) )
            println("   Next point at α = $(real(X₀[1])) ")
        
            X_grid[i] = X₀

            Dαf = DαF!( copy(Dαf₀) , component(X_grid[i],2:29)  , M₁, M₂, real(X_grid[i][1]) )
            Dₓf = DF!( copy(Dₓf₀), component(X_grid[i],2:29)  , κ, star, M₁, M₂, real(X_grid[i][1]))
            DXf = LinearOperator( space_X , space_f ,  [ Dαf[:,:]  Dₓf[:,:]] )
            
            Q, = qr(transpose(conj(Complex.(DXf[:,:]))));
            η  = Q[:,end]*conj(Q[1,end])/abs(Q[1,end])
            η[2:end] = x_remove_complex!(Sequence(space_x,η[2:end]))[:]
            η_grid[i] = Sequence(space_X,η )

            if real(η_grid[i-1] ⋅ η_grid[i]) < 0 # enforce direction
                η_grid[i] .*= -1
            end
        end
        return X_grid, η_grid, P2
    end
    X_grid1, η_grid1 , Q = PAC!( X, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂ ,N, max_arclength,direct,str_color)


    #Pseudo-Arclength Continuation 
    # Second segment
    println("\n   Second segment ")
    X2 = X_grid1[end]
    str_color = "blue"
    X_grid2, η_grid2 , Q2 = PAC!( X2, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂ ,N, max_arclength,direct,str_color)
    
    # Third segment
    println("\n   Third segment ")
    X3 = X_grid2[end]
    str_color = "blue"
    X_grid3, η_grid3 , Q3 = PAC!( X3, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂ ,N, max_arclength,direct,str_color)

#---------------------------------------------------------------------
# 7. Animation
#---------------------------------------------------------------------
    printstyled("\n7. Animation\n"; color = :blue)
    # Using data found abouve to create an animation
    # functions
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

    logocolors = Colors.JULIA_LOGO_COLORS
    colors = [logocolors.blue, logocolors.red, logocolors.green, logocolors.purple];
    Number_of_branch = 3
    for j = 1:Number_of_branch
        if j == 1
            X_grid = X_grid1
        elseif j == 2
            X_grid = X_grid2
        else
            X_grid = X_grid3
        end
        X_fft = [reverse(X_grid) ; X_grid[begin+1:end-1]]
        X_cheb = Sequence( space(X_grid[1]) , map(X -> interval.(X), grid2cheb(copy(X_fft), N)))
        α_cheb = component(X_cheb,1)[1]
        L_cheb = component(X_cheb,15)[1]
        local t = range(-1,1,100)
        line_α = zeros( 100 )
        line_L = zeros( 100 )
        for ℓ = 1:100
            line_α[ℓ] = real(mid(α_cheb(t[ℓ])))
            line_L[ℓ] = real(mid(L_cheb(t[ℓ])))
        end
        global plot_biff = plot!([line_α],[line_L], color = colors[( j + 3   ) % 4 + 1] , legend = false,markersize = 1,size=(500,500) , linewidth=3, xlabel = "α",ylabel = "L" )

        for i = 1:length(X_grid)
            local X = X_grid[i]
            local α = real(X[1])
            local x = component(X,2:29)
            local space_x = space(x)
            local γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
            scatter!([α],[L], color = colors[( j + 3   ) % 4 + 1], legend = false,markersize = 3,size=(500,500) ,markerstrokecolor =  colors[( j + 3   ) % 4 + 1], xlabel = "α",ylabel = "L")

        end
    end
    All_plots = []
    plots_biff_all = []
    for j = 1:Number_of_branch
        if j == 1
            X_grid = X_grid1
        elseif j == 2
            X_grid = X_grid2
        else
            X_grid = X_grid3
        end
        X_fft = [reverse(X_grid) ; X_grid[begin+1:end-1]]
        
        for i = 3:length(X_grid)-2
            local X = X_grid[i]
            local α = real(X[1])
            local x = component(X,2:29)
            local space_x = space(x)
            local γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
            p_temp =  plot(plot_biff) 
            scatter!([α],[L], color = "black", legend = false,markersize = 7,size=(750,750) , markerstrokecolor = "white")
            global plots_biff_all = [ plots_biff_all ; p_temp]
            global All_plots = [ All_plots ; plot_hetero_orbit2(γ,a,p,w)]
        end
    end
    titleee = strf₁*str_system( M₁ )*"\n"*strf₂*str_system( M₂ )
    anim = @animate for j = 1: length(plots_biff_all)
        plot(plots_biff_all[j] , All_plots[j] , layout=(1,2),size=(1500,750) ,plot_title = titleee ,foreground_color = "white" , 
        background_color = convert(RGB, RGB(0.18, 0.2,0.23 )) , margin=10Plots.mm )
    end
    gif(anim, "anim_GS_fps24.gif", fps = 18)

