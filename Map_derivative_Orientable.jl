# ====================================================
# Derivative of the maps
# ====================================================

function DF!(DF, x, κ, star, M₁, M₂, α)
# Components
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)

    DF[:,:] .= 0

# DΓ
    DᵧΓ = zero(component(DF,1:4,1:4)) 
    DωΓ = Sequence(codomain(component(DF,1:4,23)),component(DF,1:4,23)[:,:][:]) 
    DᵧΓ,DωΓ = DΓ!(DᵧΓ,DωΓ,γ,ω,M₁,M₂,α)  
    component(DF,1:4,1:4)[:,:] = DᵧΓ[:,:]
    component(DF,1:4,23)[:,:] = DωΓ[:]

# DV
    DᵥV = zero(component(DF,5:9,6:9)) 
    DλV = Sequence(codomain(component(DF,5:9,5)),component(DF,5:9,5)[:,:][:]) 
    DωV = Sequence(codomain(component(DF,5:9,23)),component(DF,5:9,23)[:,:][:]) 
    DᵧV = zero(component(DF,5:9,1:4)) 
    DλV, DᵥV,DωV,DᵧV  = DV!( DλV, DᵥV,DωV,DᵧV , ω, v, γ, M₁, M₂,α, λ, κ)
    component(DF,5:9,6:9)[:,:] = DᵥV[:,:]
    component(DF,5:9,5)[:,:] = DλV[:]
    component(DF,5:9,23)[:,:] = DωV[:]
    component(DF,5:9,1:4)[:,:] = DᵧV[:,:]

# DW
    DωW = Sequence(codomain(component(DF,10:13,23)),component(DF,10:13,23)[:,:][:]) 
    DλW = Sequence(codomain(component(DF,10:13,5)),component(DF,10:13,5)[:,:][:])
    DᵥW = zero(component(DF,10:13,6:9)) 
    DwW = zero(component(DF,10:13,10:13)) 
    DᵧW = zero(component(DF,10:13,1:4)) 
    DωW, DλW, DᵥW, DwW, DᵧW = DW!(DωW, DλW, DᵥW, DwW, DᵧW , w, γ, v, ω, λ, M₁, M₂,α)
    component(DF,10:13,23)[:,:] = DωW[:]
    component(DF,10:13,5)[:,:] = DλW[:]
    component(DF,10:13,6:9)[:,:] = DᵥW[:,:]
    component(DF,10:13,10:13)[:,:] = DwW[:,:]
    component(DF,10:13,1:4)[:,:] = DᵧW[:,:]

# DG
    DₐG = zero(component(DF,14:23,19:22)) 
    Dθ₀G = zero(component(DF,14:23,15:16))
    DLG = Sequence(codomain(component(DF,14:23,14)),component(DF,14:23,14)[:,:][:]) 
    DσG = zero(component(DF,14:23,17:18))
    DₚG = zero(component(DF,14:23,25:28))
    DwG = zero(component(DF,14:23,10:13))
    DₐG, Dθ₀G, DLG, DσG, DₚG,DwG = DG!(DₐG, Dθ₀G, DLG, DσG, DₚG,DwG, a, w, p, θ₀, L, σ, M₁, M₂, α)
    component(DF,14:23,19:22)[:,:] = DₐG[:,:]
    component(DF,14:23,15:16)[:,:] = Dθ₀G[:,:]
    component(DF,14:23,14)[:,:] = DLG[:]
    component(DF,14:23,17:18)[:,:] = DσG[:,:] 
    component(DF,14:23,25:28)[:,:] = DₚG[:,:]
    component(DF,14:23,10:13)[:,:] = DwG[:,:]

# DE
    component(DF,24,24)[:,:] = DE!(zero(component(DF,24,24)) , eigenpairs, M₁, M₂,α, star)[:,:]

# DP
    DeigenpairsP = zero(component(DF,25:28,24)) 
    DₚP = zero(component(DF,25:28,25:28)) 
    DeigenpairsP, DₚP = DP!(DeigenpairsP, DₚP,p, eigenpairs , M₁, M₂, α)
    component(DF,25:28,24)[:,:] = DeigenpairsP[:,:]
    component(DF,25:28,25:28)[:,:] = DₚP[:,:]


return DF
end

function DG!(DₐG, DθG, DLG, DσG, DₚG,DwG, a, w, p, θ, L, σ, M₁, M₂, α)
# Component
    a₁, a₂, a₃, a₄ = eachcomponent(a) 
    p₁, p₂, p₃, p₄ = eachcomponent(p) 
    w₁, w₂, w₃, w₄ = eachcomponent(w) 
    σ₁, σ₂ = eachcomponent(σ)
    θ₁, θ₂ = eachcomponent(θ)

# DₐG
    DₐG[:,:] .= 0
    for i = 1:4
        component(DₐG,i,i)[:,:] .= -2
        component(DₐG,i,i)[:,0] .= -1
    end
    component(DₐG,9,:)[:,:] .= 0 


    Nₐ_domain = order(domain(DₐG)[1:4])
    Nₐ_codomain = order(codomain(DₐG)[5:8])


    for i = 1:4
        for n = 0:Nₐ_domain[i]
            for m = 0:Nₐ_codomain[i]
                if n == m
                component(DₐG,i+4,i)[ m, n] = 2*n
                end
            end
        end
    end

    Da₂Ψ₁ = zeros( eltype(w), Chebyshev(Nₐ_domain[1]+1), Chebyshev(Nₐ_codomain[2] +1 )  )
    Da₁Ψ₂ = zeros( eltype(w), Chebyshev(Nₐ_domain[2]+1), Chebyshev(Nₐ_codomain[1] +1 ) )
    Da₃Ψ₂ = zeros( eltype(w), Chebyshev(Nₐ_domain[2]+1), Chebyshev(Nₐ_codomain[3] +1 ) )
    Da₄Ψ₃ = zeros( eltype(w), Chebyshev(Nₐ_domain[3]+1), Chebyshev(Nₐ_codomain[4] +1 ) )
    Da₁Ψ₄ = zeros( eltype(w), Chebyshev(Nₐ_domain[4]+1), Chebyshev(Nₐ_codomain[1] +1 ) )
    Da₃Ψ₄ = zeros( eltype(w), Chebyshev(Nₐ_domain[4]+1), Chebyshev(Nₐ_codomain[3] +1 ) )

    for i in axes(Da₂Ψ₁[:,:],1)
        for j in axes(Da₂Ψ₁[:,:],2)
            if i == j
                Da₂Ψ₁[i-1,j-1] = L/2
            end
        end
    end
    for i in axes(Da₄Ψ₃[:,:],1)
        for j in axes(Da₄Ψ₃[:,:],2)
            if i == j
                Da₄Ψ₃[i-1,j-1] = L/2
            end
        end
    end

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    if i >= 2
                    Da₁Ψ₂[:,:] = Da₁Ψ₂[:,:] + α^(ℓ-1)*L/2*M₁[i,j,ℓ] * project(  Multiplication((i-1)*a₁^(i-2)*a₃^(j-1))   ,Chebyshev(Nₐ_domain[2]+1),Chebyshev(Nₐ_codomain[1] +1 ))[:,:]
                    end
                    if j>= 2
                    Da₃Ψ₂[:,:] = Da₃Ψ₂[:,:] + α^(ℓ-1)*L/2*M₁[i,j,ℓ] * project(  Multiplication((j-1)*a₁^(i-1)*a₃^(j-2))   ,Chebyshev(Nₐ_domain[2]+1),Chebyshev(Nₐ_codomain[3] +1 ))[:,:]
                    end
                end
            end
        end
    end
    for i in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    if i >= 2
                    Da₁Ψ₄[:,:] = Da₁Ψ₄[:,:] + α^(ℓ-1)*L/2*M₂[i,j,ℓ] * project(  Multiplication((i-1)*a₁^(i-2)*a₃^(j-1))   ,Chebyshev(Nₐ_domain[4]+1),Chebyshev(Nₐ_codomain[1] +1 ))[:,:]
                    end
                    if j>= 2
                    Da₃Ψ₄[:,:] = Da₃Ψ₄[:,:] + α^(ℓ-1)*L/2*M₂[i,j,ℓ] * project(  Multiplication((j-1)*a₁^(i-1)*a₃^(j-2))   ,Chebyshev(Nₐ_domain[4]+1),Chebyshev(Nₐ_codomain[3] +1 ))[:,:]
                    end
                end
            end
        end
    end

    T = zeros( codomain(Da₂Ψ₁), codomain(Da₂Ψ₁) )
    for k = 0:order(codomain(Da₂Ψ₁))
        for j = 1:order(codomain(Da₂Ψ₁))
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end
    component(DₐG, 5,2 )[:,:] = project(T*Da₂Ψ₁, domain(component(DₐG, 5,2 )) ,codomain(component(DₐG, 5,2 )))[:,:]

    T = zeros( codomain(Da₄Ψ₃), codomain(Da₄Ψ₃) )
    for k = 0:order(codomain(Da₄Ψ₃))
        for j = 1:order(codomain(Da₄Ψ₃))
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end
    component(DₐG, 7,4 )[:,:] = project(T*Da₄Ψ₃, domain(component(DₐG, 7,4 )) ,codomain(component(DₐG, 7,4 )))[:,:]

    T = zeros( domain(Da₁Ψ₂), domain(Da₁Ψ₂) )
    for k = 0:order(domain(Da₁Ψ₂))
        for j = 1:order(domain(Da₁Ψ₂))
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end
    component(DₐG, 6,1 )[:,:] = project((T*Da₁Ψ₂), domain(component(DₐG, 6,1 )) ,codomain(component(DₐG, 6,1 )))[:,:]
    component(DₐG, 8,1 )[:,:] = project((T*Da₁Ψ₄), domain(component(DₐG, 8,1 )) ,codomain(component(DₐG, 8,1 )))[:,:]

    T = zeros( codomain(Da₃Ψ₂), codomain(Da₃Ψ₂) )
    for k = 0:order(codomain(Da₃Ψ₂))
        for j = 1:order(codomain(Da₃Ψ₂))
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end
    component(DₐG, 6,3 )[:,:] = project((T*Da₃Ψ₂), domain(component(DₐG, 6,3 )) ,codomain(component(DₐG, 6,3 )))[:,:]
    component(DₐG, 8,3 )[:,:] = project((T*Da₃Ψ₄), domain(component(DₐG, 8,3 )) ,codomain(component(DₐG, 8,3 )))[:,:]
    
    for i = 1:4
        IC = 2*ones(1,Nₐ_domain[i]+1) 
        IC[1] = 1
        IC = -IC.*((-1).^transpose(0:Nₐ_domain[i]))
        component(DₐG,i+4,i )[0,0:Nₐ_domain[i] ] = IC
    end

# Dθ₀G
        

    for i = 1:4
         N_w = order(w)[i] 
         wᵢ = component(w,i)
         n = Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
         m = Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
         component(DθG,i+4,1)[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
         component(DθG,i+4,2)[0,1] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
    end
    component(DθG,10,1)[1,1] = -2*θ₁[1]
    component(DθG,10,2)[1,1] = -2*θ₂[1]

# DLG
    DLG[:].=0
  
    Ψ₁ = 1/2*project( a₂ ,Chebyshev(order(component(DLG,5))+1) )
    Ψ₂ = Sequence(Chebyshev(order(component(DLG,6))+1),zeros(order(component(DLG,6))+2))
    Ψ₃ = 1/2*project( a₄ ,Chebyshev(order(component(DLG,7))+1) )
    Ψ₄ = Sequence(Chebyshev(order(component(DLG,8))+1),zeros(order(component(DLG,8))+2))

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ  in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    Ψ₂ = Ψ₂ + 1/2*M₁[i,j,ℓ] * α^(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₂))
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ  in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    Ψ₄ = Ψ₄ + 1/2*M₂[i,j,ℓ] * α^(ℓ-1)*project(a₁^(i-1)*a₃^(j-1),space(Ψ₄))
                end
            end 
        end
    end

    component(DLG,5)[1:end] = component(DLG,5)[1:end] + Ψ₁[2:end]  - Ψ₁[0:end-2] 
    component(DLG,6)[1:end] = component(DLG,6)[1:end] + Ψ₂[2:end]  - Ψ₂[0:end-2] 
    component(DLG,7)[1:end] = component(DLG,7)[1:end] + Ψ₃[2:end]  - Ψ₃[0:end-2] 
    component(DLG,8)[1:end] = component(DLG,8)[1:end] + Ψ₄[2:end]  - Ψ₄[0:end-2] 

    component(DLG,9)[1] =0

# DσG
    DσG[:,:] .= 0
    component( DσG ,9,1)[1,1] .= -2*σ₁
    component( DσG ,9,2)[1,1] .= -2*σ₂
    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    for i = 1:4
        component( DσG ,i,1)[1,1] = (D₁₀*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
        component( DσG ,i,2)[1,1] = 1im*(D₁₀*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
    end

# DₚG   
    DₚG[:,:] .= 0
    for i = 1:4
        Nₚ = order(domain(DₚG)[i])
        n = reshape(repeat(0:Nₚ[1],1,Nₚ[2]+1), (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:]
        m = reshape(repeat(transpose(0:Nₚ[2]),Nₚ[1]+1,1) , (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:]
        component(DₚG,i,i)[:,:] = transpose((σ₁[1] + 1im*σ₂[1] ).^n.*(σ₁[1] - 1im*σ₂[1]).^m)
    end

# DwG
    DwG[:,:] .= 0
    for i = 1:4
        N_w = order(domain(DwG))[i]    
        wᵢ = project(component(w,i) , domain(DwG)[i])
        n = Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        m = Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        component(DwG,i+4,i)[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)])[:]
        component(DwG,i+4,i)[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)])[:]
        component(DwG,i+4,i)[0,(0,:)] = ((-0.95).^m[(0,:)])[:]
    end

    return  DₐG, DθG, DLG, DσG, DₚG, DwG
end

function DP!(DeigenpairsP, DₚP,p, eigenpairs , M₁, M₂, α)
# Component
    p₁, p₂, p₃, p₄  = eachcomponent(p)
    λ₁ = component(eigenpairs,1)[1]
    ξ₁ = component(eigenpairs,2:5)[:] 
    λ₂ = component(eigenpairs,6)[1]
    ξ₂ = component(eigenpairs,7:10)[:]

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)
    M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
    M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])

# Dλ₁P
    Dλ₁P₁ = project(M₁₀*(D₁₀*p₁), codomain(DeigenpairsP)[1] )
    Dλ₁P₁[(0,0)] = 0
    Dλ₁P₁[(1,0)] = 0
    Dλ₁P₁[(0,1)] = 0
    Dλ₁P₂ = project(M₁₀*(D₁₀*p₂), codomain(DeigenpairsP)[2] )
    Dλ₁P₂[(0,0)] = 0
    Dλ₁P₂[(1,0)] = 0
    Dλ₁P₂[(0,1)] = 0
    Dλ₁P₃ = project(M₁₀*(D₁₀*p₃), codomain(DeigenpairsP)[3] )
    Dλ₁P₃[(0,0)] = 0
    Dλ₁P₃[(1,0)] = 0
    Dλ₁P₃[(0,1)] = 0
    Dλ₁P₄ = project(M₁₀*(D₁₀*p₄), codomain(DeigenpairsP)[4] )
    Dλ₁P₄[(0,0)] = 0
    Dλ₁P₄[(1,0)] = 0
    Dλ₁P₄[(0,1)] = 0

# Dλ₂P
    Dλ₂P₁ = project(M₀₁*(D₀₁*p₁), codomain(DeigenpairsP)[1] )
    Dλ₂P₁[(0,0)] = 0
    Dλ₂P₁[(1,0)] = 0
    Dλ₂P₁[(0,1)] = 0
    Dλ₂P₂ = project(M₀₁*(D₀₁*p₂), codomain(DeigenpairsP)[2] )
    Dλ₂P₂[(0,0)] = 0
    Dλ₂P₂[(1,0)] = 0
    Dλ₂P₂[(0,1)] = 0
    Dλ₂P₃ = project(M₀₁*(D₀₁*p₃), codomain(DeigenpairsP)[3] )
    Dλ₂P₃[(0,0)] = 0
    Dλ₂P₃[(1,0)] = 0
    Dλ₂P₃[(0,1)] = 0
    Dλ₂P₄ = project(M₀₁*(D₀₁*p₄), codomain(DeigenpairsP)[4] )
    Dλ₂P₄[(0,0)] = 0
    Dλ₂P₄[(1,0)] = 0
    Dλ₂P₄[(0,1)] = 0

    component(DeigenpairsP,:,1)[:,:] = [ Dλ₁P₁[:] ; Dλ₁P₂[:] ; Dλ₁P₃[:] ; Dλ₁P₄[:] ] 
    component(DeigenpairsP,:,6)[:,:] = [ Dλ₂P₁[:] ; Dλ₂P₂[:] ; Dλ₂P₃[:] ; Dλ₂P₄[:] ] 

# Dξ₁P
    component(DeigenpairsP,1:4,2:5)[:,:] .= 0 
    component(DeigenpairsP,1,2)[(1,0),:] .= -1 
    component(DeigenpairsP,2,3)[(1,0),:] .= -1 
    component(DeigenpairsP,3,4)[(1,0),:] .= -1 
    component(DeigenpairsP,4,5)[(1,0),:] .= -1 

# Dξ₂P
    component(DeigenpairsP,1:4,7:10)[:,:] .= 0 
    component(DeigenpairsP,1,7)[(0,1),:] .= -1 
    component(DeigenpairsP,2,8)[(0,1),:] .= -1 
    component(DeigenpairsP,3,9)[(0,1),:].= -1 
    component(DeigenpairsP,4,10)[(0,1),:] .= -1 

# DpP
    component(DₚP,1,1)[:,:] .=0 
    component(DₚP,1,1)[(:,1:order( component(DₚP,1,1) )[2][2]),(:,:)] = component(DₚP,1,1)[(:,1:order( component(DₚP,1,1) )[2][2]),(:,:)] + λ₂*project(D₀₁,domain(component(DₚP,1,1)),codomain(component(DₚP,1,1))[1] ⊗ Taylor( order(component(DₚP,1,1))[2][2]-1  )  )[:,:]
    component(DₚP,1,1)[(1:order( component(DₚP,1,1) )[2][1],:),(:,:)] = component(DₚP,1,1)[(1:order( component(DₚP,1,1) )[2][1],:),(:,:)] + λ₁*project(D₁₀,domain(component(DₚP,1,1)),Taylor( order(component(DₚP,1,1))[2][1]-1  ) ⊗ codomain(component(DₚP,1,1))[2]  )[:,:]
    component(DₚP,1,1)[(0,0),(:,:)] .= 0
    component(DₚP,1,1)[(0,0),(0,0)] = 1
    component(DₚP,1,1)[(1,0),(:,:)] .= 0
    component(DₚP,1,1)[(1,0),(1,0)] = 1
    component(DₚP,1,1)[(0,1),(:,:)] .= 0
    component(DₚP,1,1)[(0,1),(0,1)] = 1
    component(DₚP,1,2)[:,:] .= -1.0*Matrix(I,size(component(DₚP,1,2)[:,:]))
    component(DₚP,1,2)[(0,0),(:,:)] .= 0
    component(DₚP,1,2)[(1,0),(:,:)] .= 0
    component(DₚP,1,2)[(0,1),(:,:)] .= 0
    component(DₚP,1,3)[:,:] .=0 
    component(DₚP,1,4)[:,:] .=0 

    component(DₚP,2,1)[:,:] .=0 
    component(DₚP,2,2)[:,:] .=0 
    component(DₚP,2,2)[(:,1:order( component(DₚP,2,2) )[1][2]),(:,:)] = component(DₚP,2,2)[(:,1:order( component(DₚP,2,2) )[2][2]),(:,:)] + λ₂*project(D₀₁,domain(component(DₚP,2,2)),codomain(component(DₚP,2,2))[1] ⊗ Taylor( order(component(DₚP,2,2))[2][2]-1  )  )[:,:]
    component(DₚP,2,2)[(1:order( component(DₚP,2,2) )[1][1],:),(:,:)] = component(DₚP,2,2)[(1:order( component(DₚP,2,2) )[2][1],:),(:,:)] + λ₁*project(D₁₀,domain(component(DₚP,2,2)),Taylor( order(component(DₚP,2,2))[2][1]-1  ) ⊗ codomain(component(DₚP,2,2))[2]  )[:,:]
    component(DₚP,2,2)[(0,0),(:,:)] .= 0
    component(DₚP,2,2)[(0,0),(0,0)] = 1
    component(DₚP,2,2)[(1,0),(:,:)] .= 0
    component(DₚP,2,2)[(1,0),(1,0)] = 1
    component(DₚP,2,2)[(0,1),(:,:)] .= 0
    component(DₚP,2,2)[(0,1),(0,1)] = 1
    component(DₚP,2,3)[:,:] .=0 
    component(DₚP,2,4)[:,:] .=0 

    component(DₚP,3,1)[:,:] .=0 
    component(DₚP,3,2)[:,:] .=0 
    component(DₚP,3,3)[:,:] .=0 
    component(DₚP,3,3)[(:,1:order( component(DₚP,3,3) )[2][2]),(:,:)] = component(DₚP,3,3)[(:,1:order( component(DₚP,3,3) )[2][2]),(:,:)] + λ₂*project(D₀₁,domain(component(DₚP,3,3)),codomain(component(DₚP,3,3))[1] ⊗ Taylor( order(component(DₚP,3,3))[2][2]-1  )  )[:,:]
    component(DₚP,3,3)[(1:order( component(DₚP,3,3) )[2][1],:),(:,:)] = component(DₚP,3,3)[(1:order( component(DₚP,3,3) )[2][1],:),(:,:)] + λ₁*project(D₁₀,domain(component(DₚP,3,3)),Taylor( order(component(DₚP,3,3))[2][1]-1  ) ⊗ codomain(component(DₚP,3,3))[2]  )[:,:]
    component(DₚP,3,3)[(0,0),(:,:)] .= 0
    component(DₚP,3,3)[(0,0),(0,0)] = 1
    component(DₚP,3,3)[(1,0),(:,:)] .= 0
    component(DₚP,3,3)[(1,0),(1,0)] = 1
    component(DₚP,3,3)[(0,1),(:,:)] .= 0
    component(DₚP,3,3)[(0,1),(0,1)] = 1
    component(DₚP,3,4)[:,:] .= -1.0*Matrix(I,size(component(DₚP,3,4)[:,:]))
    component(DₚP,3,4)[(0,0),(:,:)] .= 0
    component(DₚP,3,4)[(1,0),(:,:)] .= 0
    component(DₚP,3,4)[(0,1),(:,:)] .= 0

    component(DₚP,4,1)[:,:] .=0 
    component(DₚP,4,2)[:,:] .=0 
    component(DₚP,4,3)[:,:] .=0 
    component(DₚP,4,4)[:,:] .=0 
    component(DₚP,4,4)[(:,1:order( component(DₚP,4,4) )[2][2]),(:,:)] = component(DₚP,4,4)[(:,1:order( component(DₚP,4,4) )[2][2]),(:,:)] + λ₂*project(D₀₁,domain(component(DₚP,4,4)),codomain(component(DₚP,4,4))[1] ⊗ Taylor( order(component(DₚP,4,4))[2][2]-1  )  )[:,:]
    component(DₚP,4,4)[(1:order( component(DₚP,4,4) )[2][1],:),(:,:)] = component(DₚP,4,4)[(1:order( component(DₚP,4,4) )[2][1],:),(:,:)] + λ₁*project(D₁₀,domain(component(DₚP,4,4)),Taylor( order(component(DₚP,4,4))[2][1]-1  ) ⊗ codomain(component(DₚP,4,4))[2]  )[:,:]
    component(DₚP,4,4)[(0,0),(:,:)] .= 0
    component(DₚP,4,4)[(0,0),(0,0)] = 1
    component(DₚP,4,4)[(1,0),(:,:)] .= 0
    component(DₚP,4,4)[(1,0),(1,0)] = 1
    component(DₚP,4,4)[(0,1),(:,:)] .= 0
    component(DₚP,4,4)[(0,1),(0,1)] = 1

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0 
                    if i >= 2
                        component(DₚP,2,1)[:,:] = component(DₚP,2,1)[:,:] - α^(ℓ-1)*M₁[i,j,ℓ]*(i-1)*project( Multiplication(p₁^(i-2)*p₃^(j-1)), domain(component(DₚP,2,1)),codomain(component(DₚP,2,1)))[:,:]
                    end
                    if j >= 2
                        component(DₚP,2,3)[:,:] = component(DₚP,2,3)[:,:] - α^(ℓ-1)*M₁[i,j,ℓ]*(j-1)*project( Multiplication(p₁^(i-1)*p₃^(j-2)), domain(component(DₚP,2,3)),codomain(component(DₚP,2,3)))[:,:]
                    end
                end
            end
        end
    end
    component(DₚP,2,1)[(0,0),(:,:)] .= 0
    component(DₚP,2,1)[(1,0),(:,:)] .= 0
    component(DₚP,2,1)[(0,1),(:,:)] .= 0

    component(DₚP,2,3)[(0,0),(:,:)] .= 0
    component(DₚP,2,3)[(1,0),(:,:)] .= 0
    component(DₚP,2,3)[(0,1),(:,:)] .= 0

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    if i >= 2
                        component(DₚP,4,1)[:,:] = component(DₚP,4,1)[:,:] - α^(ℓ-1)*M₂[i,j,ℓ]*(i-1)*project( Multiplication(p₁^(i-2)*p₃^(j-1)), domain(component(DₚP,4,1)),codomain(component(DₚP,4,1)))[:,:]
                    end
                    if j >= 2
                        component(DₚP,4,3)[:,:] = component(DₚP,4,3)[:,:] - α^(ℓ-1)*M₂[i,j,ℓ]*(j-1)*project( Multiplication(p₁^(i-1)*p₃^(j-2)), domain(component(DₚP,4,3)),codomain(component(DₚP,4,3)))[:,:]
                    end
                end
            end
        end
    end
    component(DₚP,4,1)[(0,0),(:,:)] .= 0
    component(DₚP,4,1)[(1,0),(:,:)] .= 0
    component(DₚP,4,1)[(0,1),(:,:)] .= 0

    component(DₚP,4,3)[(0,0),(:,:)] .= 0
    component(DₚP,4,3)[(1,0),(:,:)] .= 0
    component(DₚP,4,3)[(0,1),(:,:)] .= 0
    return DeigenpairsP, DₚP

end

function DW!(DωW, DλW, DᵥW, DwW, DᵧW , w, γ, v, ω, λ, M₁, M₂,α)
# Component
    w₁, w₂, w₃, w₄ = eachcomponent(w)
    γ_ext =  γ_extened_orientable!(γ)
    γ₁, γ₂, γ₃, γ₄ = eachcomponent(γ_ext)
    v₁, v₂, v₃, v₄ = eachcomponent(v)

    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)

# DωW
    component(DωW,1)[:] = -project(w₂, space(component(DωW,1)))[:]
    component(DωW,2)[:] .= 0
    component(DωW,3)[:] = -project(w₄, space(component(DωW,1)))[:]
    component(DωW,4)[:] .= 0
    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    component(DωW,2)[:] = component(DωW,2)[:] - α^(ℓ-1)*M₁[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(component(DωW,2)))[:]
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    component(DωW,4)[:] = component(DωW,4)[:] - α^(ℓ-1)*M₂[i,j,ℓ]*project(w₁^(i-1)*w₃^(j-1), space(component(DωW,4)))[:]
                end
            end
        end
    end

    component(DωW,1)[(:,0)] .= 0
    component(DωW,1)[(:,1)] .= 0
    component(DωW,2)[(:,0)] .= 0
    component(DωW,2)[(:,1)] .= 0
    component(DωW,3)[(:,0)] .= 0
    component(DωW,3)[(:,1)] .= 0
    component(DωW,4)[(:,0)] .= 0
    component(DωW,4)[(:,1)] .= 0
  
# DλW
    M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])

    component(DλW,1)[:] = project(M*(D₀₁*w₁),space(DλW)[1])[:]
    component(DλW,2)[:] = project(M*(D₀₁*w₂),space(DλW)[2])[:]
    component(DλW,3)[:] = project(M*(D₀₁*w₃),space(DλW)[3])[:]
    component(DλW,4)[:] = project(M*(D₀₁*w₄),space(DλW)[4])[:]

    component(DλW,1)[(:,0)] .= 0
    component(DλW,1)[(:,1)] .= 0
    component(DλW,2)[(:,0)] .= 0
    component(DλW,2)[(:,1)] .= 0
    component(DλW,3)[(:,0)] .= 0
    component(DλW,3)[(:,1)] .= 0
    component(DλW,4)[(:,0)] .= 0
    component(DλW,4)[(:,1)] .= 0

# DwW
    component(DwW,1,1)[:,:] = project(D₁₀,domain(component(DwW,1,1)),codomain(component(DwW,1,1)) )[:,:] 
    component(DwW,1,1)[(:,1:order( component(DwW,1,1) )[2][2]),(:,:)] = component(DwW,1,1)[(:,1:order( component(DwW,1,1) )[2][2]),(:,:)] + λ*project(D₀₁,domain(component(DwW,1,1)),codomain(component(DwW,1,1))[1] ⊗ Taylor( order(component(DwW,1,1))[2][2]-1  )  )[:,:]
   
    for i = -order(codomain(component(DwW,1,1)))[1]:order(codomain(component(DwW,1,1)))[1]
        for j = 0:1
            for ℓ₁ = -order(domain(component(DwW,1,1)))[1]:order(domain(component(DwW,1,1)))[1]
                for ℓ₂ = 0:order(domain(component(DwW,1,1)))[2]
                    if i == ℓ₁ && j == ℓ₂
                        component(DwW,1,1)[(i,j),(ℓ₁,ℓ₂)] = 1
                    else
                        component(DwW,1,1)[(i,j),(ℓ₁,ℓ₂)] = 0
                    end
                end
            end
        end
    end
    # component(DwW,1,1)[(:,1:order( component(DwW,1,1) )[2][2]),(:,:)] = component(DwW,1,1)[(:,1:order( component(DwW,1,1) )[2][2]),(:,:)] + λ*project(D₀₁,domain(component(DwW,1,1)),codomain(component(DwW,1,1))[1] ⊗ Taylor( order(component(DwW,1,1))[2][2]-1  )  )[:,:]
    # component(DwW,1,1)[(:,0:1),(:,:)] .= 0
    #    component(DwW,1,1)[(:,0),(:,0)] = 1.0*Matrix(I,size(component(DwW,1,1)[(:,0),(:,0)]))
    #    component(DwW,1,1)[(:,1),(:,1)] = 1.0*Matrix(I,size(component(DwW,1,1)[(:,1),(:,1)]))
    component(DwW,1,2)[:,:] .= -ω*Matrix(I,size(component(DwW,1,2)[:,:] )) 
    component(DwW,1,2)[(:,0:1),(:,:)] .= 0
    component(DwW,1,3)[:,:] .=0
    component(DwW,1,4)[:,:] .=0 

    component(DwW,2,1)[:,:].=0 
    component(DwW,2,2)[:,:].= project(D₁₀,domain(component(DwW,2,2)),codomain(component(DwW,2,2)) )[:,:]
    component(DwW,2,2)[(:,1:order( component(DwW,2,2) )[2][2]),(:,:)] = component(DwW,2,2)[(:,1:order( component(DwW,2,2) )[2][2]),(:,:)] + λ*project(D₀₁,domain(component(DwW,2,2)),codomain(component(DwW,2,2))[1] ⊗ Taylor( order(component(DwW,2,2))[2][2]-1  )  )[:,:]
    for i = -order(codomain(component(DwW,2,2)))[1]:order(codomain(component(DwW,2,2)))[1]
        for j = 0:1
            for ℓ₁ = -order(domain(component(DwW,2,2)))[1]:order(domain(component(DwW,2,2)))[1]
                for ℓ₂ = 0:order(domain(component(DwW,2,2)))[2]
                    if i == ℓ₁ && j == ℓ₂
                        component(DwW,2,2)[(i,j),(ℓ₁,ℓ₂)] = 1
                    else
                        component(DwW,2,2)[(i,j),(ℓ₁,ℓ₂)] = 0
                    end
                end
            end
        end
    end
    #   component(DwW,2,2)[(:,0:1),(:,:)] .= 0
    #   component(DwW,2,2)[(:,0),(:,0)] = 1.0*Matrix(I,size(component(DwW,2,2)[(:,0),(:,0)]))
    #   component(DwW,2,2)[(:,1),(:,1)] = 1.0*Matrix(I,size(component(DwW,2,2)[(:,1),(:,1)]))
    component(DwW,2,3)[:,:].=0 
    component(DwW,2,4)[:,:].=0 

    component(DwW,3,1)[:,:].=0 
    component(DwW,3,2)[:,:].=0 
    component(DwW,3,3)[:,:].= project(D₁₀,domain(component(DwW,3,3)),codomain(component(DwW,3,3)) )[:,:]
    component(DwW,3,3)[(:,1:order( component(DwW,3,3) )[2][2]),(:,:)] = component(DwW,3,3)[(:,1:order( component(DwW,3,3) )[2][2]),(:,:)] + λ*project(D₀₁,domain(component(DwW,3,3)),codomain(component(DwW,3,3))[1] ⊗ Taylor( order(component(DwW,3,3))[2][2]-1  )  )[:,:]
    for i = -order(codomain(component(DwW,3,3)))[1]:order(codomain(component(DwW,3,3)))[1]
        for j = 0:1
            for ℓ₁ = -order(domain(component(DwW,3,3)))[1]:order(domain(component(DwW,3,3)))[1]
                for ℓ₂ = 0:order(domain(component(DwW,3,3)))[2]
                    if i == ℓ₁ && j == ℓ₂
                        component(DwW,3,3)[(i,j),(ℓ₁,ℓ₂)] = 1
                    else
                        component(DwW,3,3)[(i,j),(ℓ₁,ℓ₂)] = 0
                    end
                end
            end
        end
    end
    #   component(DwW,3,3)[(:,0:1),(:,:)] .= 0
    #   component(DwW,3,3)[(:,0),(:,0)] = 1.0*Matrix(I,size(component(DwW,3,3)[(:,0),(:,0)]))
    #   component(DwW,3,3)[(:,1),(:,1)] = 1.0*Matrix(I,size(component(DwW,3,3)[(:,1),(:,1)]))
    component(DwW,3,4)[:,:] .= -ω*Matrix(I,size(component(DwW,3,4)[:,:] )) 
    component(DwW,3,4)[(:,0:1),(:,:)] .= 0

    component(DwW,4,1)[:,:].=0 
    component(DwW,4,2)[:,:].=0 
    component(DwW,4,3)[:,:].=0 
    component(DwW,4,4)[:,:].= project(D₁₀,domain(component(DwW,4,4)),codomain(component(DwW,4,4)) )[:,:]
    component(DwW,4,4)[(:,1:order( component(DwW,4,4) )[2][2]),(:,:)] = component(DwW,4,4)[(:,1:order( component(DwW,4,4) )[2][2]),(:,:)] + λ*project(D₀₁,domain(component(DwW,4,4)),codomain(component(DwW,4,4))[1] ⊗ Taylor( order(component(DwW,4,4))[2][2]-1  )  )[:,:]
    for i = -order(codomain(component(DwW,4,4)))[1]:order(codomain(component(DwW,4,4)))[1]
        for j = 0:1
            for ℓ₁ = -order(domain(component(DwW,4,4)))[1]:order(domain(component(DwW,4,4)))[1]
                for ℓ₂ = 0:order(domain(component(DwW,4,4)))[2]
                    if i == ℓ₁ && j == ℓ₂
                        component(DwW,4,4)[(i,j),(ℓ₁,ℓ₂)] = 1
                    else
                        component(DwW,4,4)[(i,j),(ℓ₁,ℓ₂)] = 0
                    end
                end
            end
        end
    end
    #   component(DwW,4,4)[(:,0:1),(:,:)] .= 0
    #    component(DwW,4,4)[(:,0),(:,0)] = 1.0*Matrix(I,size(component(DwW,4,4)[(:,0),(:,0)]))
    #    component(DwW,4,4)[(:,1),(:,1)] = 1.0*Matrix(I,size(component(DwW,4,4)[(:,1),(:,1)]))

    for i = 1 :4
        for k = -order(codomain(component(DwW,i,i))[1]):order(codomain(component(DwW,i,i))[1])
            for ℓ = -order(domain(component(DwW,i,i))[1]):order(domain(component(DwW,i,i))[1])
                if ℓ == k
                    component(DwW,i,i)[(ℓ,0),(k,0)] = convert(eltype(w),1.0)
                    component(DwW,i,i)[(ℓ,1),(k,1)] = convert(eltype(w),1.0)
                end
            end
        end
    end

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0 
                    if i >= 2
                        component(DwW,2,1)[:,:] = component(DwW,2,1)[:,:] - α^(ℓ-1)*ω*M₁[i,j,ℓ]*(i-1)*project( Multiplication(w₁^(i-2)*w₃^(j-1)), domain(component(DwW,2,1)) ,codomain(component(DwW,2,1)))[:,:]
                    end
                    if j >= 2
                        component(DwW,2,3)[:,:] = component(DwW,2,3)[:,:] - α^(ℓ-1)*ω*M₁[i,j,ℓ]*(j-1)*project( Multiplication(w₁^(i-1)*w₃^(j-2)), domain(component(DwW,2,3)) ,codomain(component(DwW,2,3)))[:,:]
                    end
                end
            end
        end
    end
    component(DwW,2,1)[(:,0:1),(:,:)] .= 0
    component(DwW,2,3)[(:,0:1),(:,:)] .= 0


    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0
                    if i >= 2
                        component(DwW,4,1)[:,:] = component(DwW,4,1)[:,:] - α^(ℓ-1)*ω*M₂[i,j,ℓ]*(i-1)*project( Multiplication(w₁^(i-2)*w₃^(j-1)), domain(component(DwW,4,1)) ,codomain(component(DwW,4,1)))[:,:]
                    end
                    if j >= 2
                        component(DwW,4,3)[:,:] = component(DwW,4,3)[:,:] - α^(ℓ-1)*ω*M₂[i,j,ℓ]*(j-1)*project( Multiplication(w₁^(i-1)*w₃^(j-2)), domain(component(DwW,4,3)) ,codomain(component(DwW,4,3)))[:,:]
                    end
                end
            end
        end
    end
    component(DwW,4,1)[(:,0:1),(:,:)] .= 0
    component(DwW,4,3)[(:,0:1),(:,:)] .= 0

# DᵥW
    DᵥW[:,:] .= 0 
    for i = 1 :4
        for k = -order(domain(component(DᵥW,i,i))):order(domain(component(DᵥW,i,i)))
            for ℓ = -order(codomain(component(DᵥW,i,i))[1]):order(codomain(component(DᵥW,i,i))[1])
                if ℓ == k
                    component(DᵥW,i,i)[(ℓ,1),k] = convert(eltype(w),-1.0)
                end
            end
        end
    end
    # component(DᵥW,1,1)[(:,1),:] = (-convert.(eltype(w),1.0*Matrix(I,2*order(codomain(component(DᵥW,1,1))[1])+1,2*order(domain(component(DᵥW,1,1)))+1 )))
    #component(DᵥW,2,2)[(:,1),:] = (-convert.(eltype(w),1.0*Matrix(I,2*order(codomain(component(DᵥW,2,2))[1])+1,2*order(domain(component(DᵥW,2,2)))+1 )))
    #component(DᵥW,3,3)[(:,1),:] = (-convert.(eltype(w),1.0*Matrix(I,2*order(codomain(component(DᵥW,3,3))[1])+1,2*order(domain(component(DᵥW,3,3)))+1 )))
    #component(DᵥW,4,4)[(:,1),:] = (-convert.(eltype(w),1.0*Matrix(I,2*order(codomain(component(DᵥW,4,4))[1])+1,2*order(domain(component(DᵥW,4,4)))+1 )))

# DᵧW
    DᵧW[:,:] .=0
    for i = 1:4
        Nw = order(codomain(DᵧW)[i])[1]
        Nt = order(domain(DᵧW)[i])[1]
        if isodd(i)
            D_even = zeros(eltype(w),2*Nw+1 , Nt+1)
            for j = 0:min(Nw,Nt)
                if j == 0
                    D_even[1+Nw, j+1] = 1
                else
                    D_even[1+Nw+j, j+1] = 1
                    D_even[1+Nw-j, j+1] = 1
                end
            end
            component( DᵧW ,i,i)[(:,0),:] = -D_even
        else
            D_odd= zeros(eltype(w),2*Nw +1,  Nt )
            for j = 1:min(Nw,Nt)
                D_odd[Nw+j+1, j] =-1im
                D_odd[Nw-j+1, j] = 1im
            end
            component( DᵧW,i,i)[(:,0),:] = -D_odd
        end

    end
    return  DωW, DλW, DᵥW, DwW, DᵧW
end

function DE!(DE, eigenpairs, M₁, M₂,α, star)
    if length(star) == 1
        star = [star; star]
    end
# DE 
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
    M = ExactReal.([0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0])
    for  ℓ in axes(M₁,3)
        M = M + [ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); α^ExactReal(ℓ-1)*M₁_temp[2,1,ℓ] ExactReal(0) α^ExactReal(ℓ-1)*M₁_temp[1,2,ℓ] ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0)]
    end
    for  ℓ in axes(M₂,3)
        M = M + [ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); ExactReal(0) ExactReal(0) ExactReal(0) ExactReal(0); α^ExactReal(ℓ-1)*M₂_temp[2,1,ℓ] ExactReal(0) α^ExactReal(ℓ-1)*M₂_temp[1,2,ℓ] ExactReal(0)]
    end

    δₖ₁ = ExactReal(0)
    δₖ₂ = ExactReal(0)
    δₖ₃ = ExactReal(0)
    δₖ₄ = ExactReal(0)

    if star[1] == 1
        δₖ₁ = ExactReal(1)
    elseif  star[1] == 2
        δₖ₂ = ExactReal(1)
    elseif  star[1] == 3
        δₖ₃ = ExactReal(1)
    elseif  star[1] == 4
        δₖ₄ = ExactReal(1)
    end

    Δₖ₁ = ExactReal(0)
    Δₖ₂ = ExactReal(0)
    Δₖ₃ = ExactReal(0)
    Δₖ₄ = ExactReal(0)
    if star[2] == 1
        Δₖ₁ = ExactReal(1)
    elseif  star[2] == 2
        Δₖ₂ = ExactReal(1)
    elseif  star[2] == 3
        Δₖ₃ = ExactReal(1)
    elseif  star[2] == 4
        Δₖ₄ = ExactReal(1)
    end

    if eltype(α) == Interval{Float64}
        DE[1,:] = interval.([ 0 δₖ₁ δₖ₂ δₖ₃ δₖ₄ 0 0 0 0 0])
        DE[2:5,:] = [ξ₁ interval.(1.0*Matrix(I,4,4))*λ₁-M interval.(zeros(4,5)) ]   
        DE[6,:] = interval.([ 0 0 0 0 0 0 Δₖ₁ Δₖ₂ Δₖ₃ Δₖ₄])
        DE[7:10,:] = [interval.(zeros(4,5)) ξ₂ interval.(1.0*Matrix(I,4,4))*λ₂-M ]  
    else
        DE[1,:] = [ 0 δₖ₁ δₖ₂ δₖ₃ δₖ₄ 0 0 0 0 0]
        DE[2:5,:] = [ξ₁ 1.0*Matrix(I,4,4)*λ₁-M zeros(4,5) ]   
        DE[6,:] = [ 0 0 0 0 0 0 Δₖ₁ Δₖ₂ Δₖ₃ Δₖ₄]
        DE[7:10,:] = [zeros(4,5) ξ₂ 1.0*Matrix(I,4,4)*λ₂-M ]  
    end
    return DE

end


function DV!( DλV, DᵥV,DωV,DᵧV , ω, v, γ, M₁, M₂,α, λ, κ)
# Parameters 
    γ_ext =  γ_extened_orientable!(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    γ₁ = component(γ_ext,1)
    γ₃ = component(γ_ext,3)
    D = Derivative(1)
    
# DωV 
    component(DωV,1)[1] = 0
    component(DωV,2)[:] = project(-v₂,space(DωV)[2])[:]
    component(DωV,3)[:] .= 0
    component(DωV,4)[:] = project(-v₄,space(DωV)[4])[:]
    component(DωV,5)[:] .= 0

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n-2 >= 0 && m-2>= 0
                        component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[3])[:]
                    elseif n-2 >= 0
                        component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ ,space(DωV)[3])[:]
                    elseif m-2 >= 0
                        component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[3])[:]
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
                        component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃  ,space(DωV)[5])[:]
                    elseif n-2 >= 0
                        component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁  ,space(DωV)[5])[:]
                    elseif m-2 >= 0
                        component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[5])[:]
                    end
                end
            end
        end
    end

# DλV
    component(DλV,1)[1]= 0.0
    component(DλV,2)[:] = project(v₁, space(DλV)[2] )[:]
    component(DλV,3)[:] = project(v₂, space(DλV)[3] )[:]
    component(DλV,4)[:] = project(v₃, space(DλV)[4] )[:]
    component(DλV,5)[:] = project(v₄, space(DλV)[5] )[:]
 
# DᵥV
    component(DᵥV,2,1)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[1],codomain(DᵥV)[2], eltype(DᵥV) )[:,:]) + project(Multiplication(λ*(v₁^0)),domain(DᵥV)[1],codomain(DᵥV)[2])[:,:] 
    component(DᵥV,2,2)[:,:] = -ω*Matrix(I, 2*order(codomain(DᵥV)[2])+1, 2*order(domain(DᵥV)[2])+1)
    component(DᵥV,2,3)[:,:] = zeros(2*order(codomain(DᵥV)[2])+1,2*order(domain(DᵥV)[3])+1)
    component(DᵥV,2,4)[:,:] = zeros(2*order(codomain(DᵥV)[2])+1,2*order(domain(DᵥV)[4])+1)

    component(DᵥV,3,1)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[1])+1)
    component(DᵥV,3,2)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[2],codomain(DᵥV)[3], eltype(DᵥV))[:,:]) + project(Multiplication(λ*(v₂^0)),domain(DᵥV)[2],codomain(DᵥV)[3])[:,:] 
    component(DᵥV,3,3)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[3])+1)
    component(DᵥV,3,4)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[4])+1)

    component(DᵥV,4,1)[:,:] = zeros(2*order(codomain(DᵥV)[4])+1,2*order(domain(DᵥV)[1])+1)
    component(DᵥV,4,2)[:,:] = zeros(2*order(codomain(DᵥV)[4])+1,2*order(domain(DᵥV)[2])+1)
    component(DᵥV,4,3)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[3],codomain(DᵥV)[4], eltype(DᵥV) )[:,:]) + project(Multiplication(λ*(v₃^0)), domain(DᵥV)[3],codomain(DᵥV)[4])[:,:] 
    component(DᵥV,4,4)[:,:] = -ω*Matrix(I, 2*order(codomain(DᵥV)[4])+1, 2*order(domain(DᵥV)[4])+1)

    component(DᵥV,5,1)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[1])+1)
    component(DᵥV,5,2)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[2])+1)
    component(DᵥV,5,3)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[3])+1)
    component(DᵥV,5,4)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[4],codomain(DᵥV)[5], eltype(DᵥV))[:,:]) + project(Multiplication(λ*(v₄^0)), domain(DᵥV)[4],codomain(DᵥV)[5])[:,:] 

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n-2 >= 0
                        component(DᵥV,3,1)[:,:] = component(DᵥV,3,1)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ) ,domain(DᵥV)[1], codomain(DᵥV)[3])[:,:]
                    end
                    if m-2 >= 0
                        component(DᵥV,3,3)[:,:] = component(DᵥV,3,3)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DᵥV)[3], codomain(DᵥV)[3])[:,:]
                    end
                end
            end
        end
    end

    for n in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    if n-2 >= 0
                        component(DᵥV,5,1)[:,:] = component(DᵥV,5,1)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ),domain(DᵥV)[1], codomain(DᵥV)[5])[:,:]
                    end
                    if m-2 >= 0
                        component(DᵥV,5,3)[:,:] = component(DᵥV,5,3)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DᵥV)[3], codomain(DᵥV)[5])[:,:]
                    end
                end
            end
        end
    end
    component(DᵥV,1,1)[:,:] .= 0
    component(DᵥV,1,2)[:,:] .= 0
    component(DᵥV,1,3)[:,:] .= 0
    component(DᵥV,1,4)[:,:] .= 0
    component(DᵥV,1,κ)[:,:] .= 1
 
    
    
# DᵧV
    DᵧV[:,:] .= 0

    codomainᵥ = codomain(DᵧV)[2:5]
    domainᵥ = domain(DᵧV)
    Nᵥ =  order(codomainᵥ)
    Nᵧ =  order(domainᵥ)

    Dγ₁V₂temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[2])
    Dγ₃V₂temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[2])
    Dγ₁V₄temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[4])
    Dγ₃V₄temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[4])

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n-3 >= 0 
                        Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  -ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*v₁ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:] 
                    end
                    if n-2 >= 0 && m-2 >= 0
                        Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₃ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:]
                        Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₁ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:]  
                    end
                    if m-3 >=0
                        Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:] 
                    end
                end
            end
        end
    end

    for n in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    if n-3 >= 0 
                        Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*v₁ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:] 
                    end
                    if n-2 >= 0 && m-2 >= 0
                        Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₃ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:]
                        Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₁ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:]  
                    end
                    if m-3 >=0
                        Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:] 
                    end
                end
            end
        end
    end


    D_even = zeros(eltype(Dγ₁V₂temp) ,2*Nᵧ[1]+1 , Nᵧ[1] +1)
    for j = 0:Nᵧ[1]
        if j == 0
            D_even[1+Nᵧ[1], j+1] = 1
        else
            D_even[1+Nᵧ[1]+j, j+1] = 1
            D_even[1+Nᵧ[1]-j, j+1] = 1
        end
    end
    component(DᵧV, 3 ,1)[:,:] = Dγ₁V₂temp[:,:]*D_even[:,:]
    component(DᵧV, 3 ,3)[:,:] = Dγ₃V₂temp[:,:]*D_even[:,:]
    
    
    D_even = zeros(eltype(Dγ₁V₂temp),2*Nᵧ[3]+1 , Nᵧ[3] +1)
    for j = 0:Nᵧ[3]
        if j == 0
            D_even[1+Nᵧ[3], j+1] = 1
        else
            D_even[1+Nᵧ[3]+j, j+1] = 1
            D_even[1+Nᵧ[3]-j, j+1] = 1
        end
    end
    component(DᵧV, 5 ,1)[:,:] = Dγ₁V₄temp[:,:]*D_even[:,:]
    component(DᵧV, 5 ,3)[:,:] = Dγ₃V₄temp[:,:]*D_even[:,:]   
    return  DλV, DᵥV, DωV, DᵧV 
end

function DΓ!(DᵧΓ,DωΓ,γ,ω,M₁,M₂,α)  
# Parameters
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
    
# DωΓ
    D = Derivative(1)
    component(DωΓ,1)[:]  = project( γ₂, space(DωΓ)[1] )[:]
    component(DωΓ,2)[:] .= 0
    component(DωΓ,3)[:]  = project( γ₄, space(DωΓ)[3] )[:]
    component(DωΓ,4)[:] .= 0
    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    component(DωΓ,2)[:] = component(DωΓ,2)[:] - α^(ℓ-1)*project(M₁[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(DωΓ)[2])[:]
                end
            end
        end
    end
    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    component(DωΓ,4)[:] = component(DωΓ,4)[:] - α^(ℓ-1)*project(M₂[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), space(DωΓ)[4])[:]
                end
            end
        end
    end
    
# DᵧΓ
    component(DᵧΓ,1,1)[:,:] = -Matrix(project(Derivative(1),domain(DᵧΓ)[1],codomain(DᵧΓ)[1], eltype(DᵧΓ) )[:,:])
    component(DᵧΓ,1,2)[:,:] = ω*Matrix( I, order(codomain(DᵧΓ)[1]), order(domain(DᵧΓ)[2]))
    component(DᵧΓ,1,3)[:,:] = zeros(order(codomain(DᵧΓ)[1]),order(domain(DᵧΓ)[3])+1)
    component(DᵧΓ,1,4)[:,:] = zeros(order(codomain(DᵧΓ)[1]),order(domain(DᵧΓ)[4]))

    component(DᵧΓ,2,1)[:,:] = zeros(order(codomain(DᵧΓ)[2])+1,order(domain(DᵧΓ)[1])+1)
    component(DᵧΓ,2,2)[:,:] = Matrix(project(Derivative(1),domain(DᵧΓ)[2],codomain(DᵧΓ)[2], eltype(DᵧΓ))[:,:])
    component(DᵧΓ,2,3)[:,:] = zeros(order(codomain(DᵧΓ)[2])+1,order(domain(DᵧΓ)[3])+1)
    component(DᵧΓ,2,4)[:,:] = zeros(order(codomain(DᵧΓ)[2])+1,order(domain(DᵧΓ)[4]))

    component(DᵧΓ,3,1)[:,:] = zeros(order(codomain(DᵧΓ)[3]),order(domain(DᵧΓ)[1])+1)
    component(DᵧΓ,3,2)[:,:] = zeros(order(codomain(DᵧΓ)[3]),order(domain(DᵧΓ)[2]))
    component(DᵧΓ,3,3)[:,:] = -Matrix(project(Derivative(1),domain(DᵧΓ)[3],codomain(DᵧΓ)[3], eltype(DᵧΓ) )[:,:])
    component(DᵧΓ,3,4)[:,:] = ω*Matrix(I, order(codomain(DᵧΓ)[3]), order(domain(DᵧΓ)[4]))

    component(DᵧΓ,4,1)[:,:] = zeros(order(codomain(DᵧΓ)[4])+1,order(domain(DᵧΓ)[1])+1)
    component(DᵧΓ,4,2)[:,:] = zeros(order(codomain(DᵧΓ)[4])+1,order(domain(DᵧΓ)[2]))
    component(DᵧΓ,4,3)[:,:] = zeros(order(codomain(DᵧΓ)[4])+1,order(domain(DᵧΓ)[3])+1)
    component(DᵧΓ,4,4)[:,:] = Matrix(project(Derivative(1),domain(DᵧΓ)[4],codomain(DᵧΓ)[4], eltype(DᵧΓ) )[:,:])

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0 
                    if n >= 2
                        component(DᵧΓ,2,1)[:,:] = component(DᵧΓ,2,1)[:,:] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(n-1)*project( Multiplication(γ₁^(n-2)*γ₃^(m-1)),domain(DᵧΓ)[1],codomain(DᵧΓ)[2])[:,:]
                    end
                    if m >= 2
                        component(DᵧΓ,2,3)[:,:] = component(DᵧΓ,2,3)[:,:] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(m-1)*project( Multiplication(γ₁^(n-1)*γ₃^(m-2)),domain(DᵧΓ)[3],codomain(DᵧΓ)[2])[:,:]
                    end
                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    if n >= 2
                        component(DᵧΓ,4,1)[:,:] = component(DᵧΓ,4,1)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(n-1)*project( Multiplication(γ₁^(n-2)*γ₃^(m-1)),domain(DᵧΓ)[1],codomain(DᵧΓ)[4])[:,:]
                    end
                    if m >= 2
                        component(DᵧΓ,4,3)[:,:] = component(DᵧΓ,4,3)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(m-1)*project( Multiplication(γ₁^(n-1)*γ₃^(m-2)),domain(DᵧΓ)[3],codomain(DᵧΓ)[4])[:,:]
                    end
                end
            end
        end
    end

    return DᵧΓ,DωΓ

end

function DαF!( DαF, x, M₁, M₂, α )
# DαF
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    w₁,w₂,w₃,w₄ = eachcomponent(w)
    p₁, p₂, p₃, p₄ = eachcomponent(p)
    a₁, a₂, a₃, a₄ = eachcomponent(a)

    γ_ext = γ_extened_orientable!(γ)
    γ₁_ext = component(γ_ext,1)
    γ₃_ext = component(γ_ext,3)


    DαΨ₂ = Sequence(Chebyshev(order( codomain(component(DαF,19)) )+1),zeros(eltype(x),order(codomain(component(DαF,19)) )+2))
    DαΨ₄ = Sequence(Chebyshev(order( codomain(component(DαF,21)) )+1),zeros(eltype(x),order(codomain(component(DαF,21)) )+2))

    #DαΨ₂ = Sequence(Chebyshev(order(component(x,19))+1),zeros(eltype(x),order(component(x,19))+2))
    #DαΨ₄ = Sequence(Chebyshev(order(component(x,21))+1),zeros(eltype(x),order(component(x,21))+2))


    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0 && ℓ > 1
                    component(DαF,2)[:,:]  = component(DαF,2)[:,:]  - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * project(γ₁^(n-1)*γ₃^(m-1),codomain(component(DαF,2 )))[:]
                    component(DαF,11)[:,:] = component(DαF,11)[:,:]- (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * project(w₁^(n-1)*w₃^(m-1),codomain(component(DαF,11)))[:]
                    component(DαF,26)[:,:] = component(DαF,26)[:,:] - (ℓ-1)*α^(ℓ-2)*M₁[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),codomain(component(DαF,26)))[:]
                    DαΨ₂ = DαΨ₂ + L/2*(ℓ-1)*α^(ℓ-2)*M₁[n,m,ℓ] * project(a₁^(n-1)*a₃^(m-1),space(DαΨ₂))
                    if n-2 >= 0 && m-2>= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * project((n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(component(DαF,7)))[:]
                    elseif n-2 >= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * project((n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ ,codomain(component(DαF,7)))[:]
                    elseif m-2 >= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * project((m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(component(DαF,7)))[:]
                    end
                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0 && ℓ > 1
                    component(DαF,4 )[:,:] = component(DαF,4 )[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * project(γ₁^(n-1)*γ₃^(m-1),codomain(component(DαF,4 )))[:]
                    component(DαF,13)[:,:] = component(DαF,13)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * project(w₁^(n-1)*w₃^(m-1),codomain(component(DαF,13)))[:]
                    component(DαF,28)[:,:] = component(DαF,28)[:,:] - (ℓ-1)*α^(ℓ-2)*M₂[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),codomain(component(DαF,28)))[:]
                    DαΨ₄ = DαΨ₄ + L/2*(ℓ-1)*α^(ℓ-2)*M₂[n,m,ℓ] * project(a₁^(n-1)*a₃^(m-1),space(DαΨ₄))
                    if n-2 >= 0 && m-2>= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * project((n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(component(DαF,9)))[:]
                    elseif n-2 >= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * project((n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ ,codomain(component(DαF,9)))[:]
                    elseif m-2 >= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * project((m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(component(DαF,9)))[:]
                    end
                end
            end
        end
    end

    component(DαF,11)[(:,0),1] .= 0
    component(DαF,11)[(:,1),1] .= 0
    component(DαF,13)[(:,0),1] .= 0
    component(DαF,13)[(:,1),1] .= 0

    component(DαF,26)[(0,0),1] = 0
    component(DαF,26)[(1,0),1] = 0
    component(DαF,26)[(0,1),1] = 0
    component(DαF,28)[(0,0),1] = 0
    component(DαF,28)[(1,0),1] = 0
    component(DαF,28)[(0,1),1] = 0

    component(DαF,19)[1:end,1] = DαΨ₂[2:end] - DαΨ₂[0:end-2] 
    component(DαF,21)[1:end,1] = DαΨ₄[2:end] - DαΨ₄[0:end-2] 
    component(DαF,19)[0,1] = 0
    component(DαF,21)[0,1] = 0


    M = zeros(eltype(x),4,4)
    for ℓ in axes(M₁,3)
        if ℓ > 1 
            if size(M₁,1) > 1
                M[2,1] = M[2,1] + M₁[2,1,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
            if size(M₁,2) > 1
                M[2,3] = M[2,3] + M₁[1,2,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
        end
    end
    for ℓ in axes(M₂,3)
        if ℓ > 1 
            if size(M₂,1) > 1
                M[4,1] = M[4,1] + M₂[2,1,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
            if size(M₂,2) > 1
                M[4,3] = M[4,3] + M₂[1,2,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
        end
    end
    component(DαF,24)[:,:] =  [ 0;- M*component(eigenpairs,2:5)[:]; 0; - M*component(eigenpairs,7:10)[:] ]

     return DαF
end

function DF_all!(DF, X, κ, star, lengthᵥ, M₁, M₂ , Xₖ_dot )
    α = real(X[1])
    x = component(X,2:29)

    component(DF,2:29,1)[:,:] =  DαF!( component(DF,2:29,1), x, M₁, M₂, α )[:,:]

    component(DF,2:29,2:29)[:,:]  =  DF!(component(DF,2:29,2:29), x, κ, star, M₁, M₂, α)[:,:]
    #component(DF,1,:)[:,:] =  Xₖ_dot[:]
    component(DF,1,:)[:,:] =  conj(project(Xₖ_dot, domain(component(DF,1,:))  ) )[:]
    return DF

end

function DαΓ!( DαF, γ , M₁, M₂, α, ω )
    # DαF
        γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
        DαF[:] .= 0
        for i in axes(M₁,1)
            for j in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[i,j,ℓ] != 0 && ℓ > 1
                        component(DαF,2)[:]  = component(DαF,2)[:]  - (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project(γ₁^(i-1)*γ₃^(j-1),space(component(DαF,2 )))[:]
                    end
                end
            end
        end
    
        for i  in axes(M₂,1)
            for j in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[i,j,ℓ] != 0 && ℓ > 1
                        component(DαF,4 )[:] = component(DαF,4 )[:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project(γ₁^(i-1)*γ₃^(j-1),space(component(DαF,4 )))[:]
                    end
                end
            end
        end
        return DαF
end
    

function DΓ_all!(DF, X, ω, M₁, M₂ , Xₖ_dot )
    α = real(component(X,1)[1])
    γ = component(X,2:5)

    DF[:,:] .= 0
    DαΓ = zeros(space(Derivative(1)*γ))
    component(DF,2:5,1)[:,:] =  DαΓ!( DαΓ, γ , M₁, M₂, α, ω )[:]
    
    Df  = LinearOperator( space(γ), space(Derivative(1)*γ) , zeros(length(DαΓ), length(γ) )  )  
    DωΓ = zeros(space(Derivative(1)*γ))
    Df, DωΓ  =  DΓ!(Df,DωΓ,γ,ω,M₁,M₂,α)
    component(DF,2:5,2:5)[:,:] = Df[:,:]
    #component(DF,1,:)[:,:] =  Xₖ_dot[:]
    component(DF,1,:)[:,:] =  conj(Xₖ_dot)[:]
    return DF

end


function DFᵢⱼ(DF, X, κ, star, lengthᵥ, M₁, M₂ , Xₖ_dot, i, j )
# Variables
    α = real(X[1])
    x = component(X,2:29)
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
    λ₁ = component( eigenpairs , 1)[1]
    λ₂ = component( eigenpairs , 6)[1]
    γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
    v₁,v₂,v₃,v₄ = eachcomponent(v)
    w₁,w₂,w₃,w₄ = eachcomponent(w)
    p₁, p₂, p₃, p₄ = eachcomponent(p)
    a₁, a₂, a₃, a₄ = eachcomponent(a)
    σ₁, σ₂ = eachcomponent(σ)
    θ₁, θ₂ = eachcomponent(θ₀)

    γ_ext = γ_extened_orientable!(γ)
    γ₁_ext = component(γ_ext,1)
    γ₃_ext = component(γ_ext,3)
# i = 1 =============================================================================
    if i == 1
        DF[:,:] =  conj(project( component( Xₖ_dot ,j) , domain(DF)  ) )[:]
# i = 2 =============================================================================
    elseif i == 2
        if j == 2
            DF[:,:]  = -Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF) )[:,:])
        elseif j == 3
            if eltype(α) ==  Interval{Float64}
                DF[:,:]  = ω*interval.(Matrix( I , order(codomain(DF)), order(domain(DF))))
            else
                DF[:,:]  = ω*Matrix( I , order(codomain(DF)), order(domain(DF)))
            end
        elseif j == 24
            DF[:,:]  = project( γ₂, codomain(DF) )[:]
        end
# i = 3 =============================================================================
    elseif i == 3
        if j == 1
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:]  = DF[:,:]  - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₁[n,m,ℓ] * project(γ₁^(n-1)*γ₃^(m-1),codomain(DF))[:]
                        end
                    end
                end
            end
        elseif j == 2
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(γ₁^(n-2)*γ₃^(m-1)),domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 3
            DF[:,:] = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF))[:,:])
        elseif j == 4 
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(γ₁^(n-1)*γ₃^(m-2)),domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end

        elseif j == 24
            DF[:,:]  .= ExactReal(0)
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*project(M₁[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), codomain(DF))[:]
                        end
                    end
                end
            end
        end
# i = 4 =============================================================================
    elseif i == 4
        if j == 4
            DF[:,:]  = -Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF) )[:,:])
        elseif j == 5
            if eltype(α) ==  Interval{Float64}
                DF[:,:]  = ω*interval.(Matrix( I , order(codomain(DF)), order(domain(DF))))
            else
                DF[:,:]  = ω*Matrix( I , order(codomain(DF)), order(domain(DF)))
            end
        elseif j == 24
            DF[:,:]  = project( γ₄, codomain(DF) )[:]
        end
# i = 5 =============================================================================
    elseif i == 5
        if j == 1
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₂[n,m,ℓ] * project(γ₁^(n-1)*γ₃^(m-1),codomain(DF))[:]
                        end
                    end
                end
            end
        
        elseif j == 2
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(γ₁^(n-2)*γ₃^(m-1)),domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 4
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(γ₁^(n-1)*γ₃^(m-2)),domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 5
            DF[:,:] = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF) )[:,:])
        elseif j == 24
            DF[:,:]  .= ExactReal(0)
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            DF[:,:]  = DF[:,:]  - α^ExactReal(ℓ-1)*project(M₂[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)), codomain(DF))[:]
                        end
                    end
                end
            end
        end
# i = 6 =============================================================================
    elseif i == 6
        if j == 6 + κ 
            DF[:,:] .= ExactReal(1)
        end
# i = 7 =============================================================================
    elseif i == 7
        if j == 6
            DF[:,:] = project(v₁, codomain(DF))[:]
        elseif j == 7
            DF[:,:]  = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF) )[:,:]) + project(Multiplication(λ*(v₁^0)),domain(DF),codomain(DF))[:,:] 
        elseif j ==8
            if eltype(α) ==  Interval{Float64}
                DF[:,:]  = -ω*interval.(Matrix(I, 2*order(codomain(DF))+1, 2*order(domain(DF))+1))
            else
                DF[:,:]  = -ω*Matrix(I, 2*order(codomain(DF))+1, 2*order(domain(DF))+1)
            end
        elseif j == 24 
            DF[:,:] = project(-v₂,codomain(DF))[:]
        end
# i = 8 =============================================================================
    elseif i == 8
        if j == 1
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 && ℓ > 1
                            if n-2 >= 0 && m-2>= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₁[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            elseif n-2 >= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₁[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ ,codomain(DF))[:]
                            elseif m-2 >= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₁[n,m,ℓ] * project(ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            end
                        end
                    end
                end
            end
        elseif j == 2
                    
            codomainᵥ = codomain(DF)
            domainᵥ = domain(DF)
            Nᵥ =  order(codomainᵥ)
            Nᵧ =  order(domainᵥ)

            DγVtemp = zeros(eltype(DF), Fourier(order(domainᵥ), interval(1.0) ),codomainᵥ)


            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if n-3 >= 0 
                                DγVtemp[:,:] = DγVtemp[:,:]  -ω*α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( ExactReal((n-1)*(n-2))*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:] 
                            end
                            if n-2 >= 0 && m-2 >= 0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication(ExactReal( (m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:]
                            end

                        end
                    end
                end
            end

            D_even = zeros(eltype(DγVtemp) ,2*Nᵧ+1 , Nᵧ +1)
            for n = 0:Nᵧ
                if n == 0
                    D_even[1+Nᵧ, n+1] = ExactReal(1)
                else
                    D_even[1+Nᵧ+n, n+1] = ExactReal(1)
                    D_even[1+Nᵧ-n, n+1] = ExactReal(1)
                end
            end
            DF[:,:] = DγVtemp[:,:]*D_even[:,:]


        elseif j == 4
            codomainᵥ = codomain(DF)
            domainᵥ = domain(DF)
            Nᵥ =  order(codomainᵥ)
            Nᵧ =  order(domainᵥ)

            DγVtemp = zeros(eltype(DF), Fourier(order(domainᵥ), interval(1.0) ),codomainᵥ)

            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if n-2 >= 0 && m-2 >= 0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:]  
                            end
                            if m-3 >=0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( ExactReal((m-1)*(m-2))*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:] 
                            end
                        end
                    end
                end
            end
            D_even = zeros(eltype(DγVtemp) ,2*Nᵧ+1 , Nᵧ +1)
            for n = 0:Nᵧ
                if n == 0
                    D_even[1+Nᵧ, n+1] = ExactReal(1)
                else
                    D_even[1+Nᵧ+n, n+1] = ExactReal(1)
                    D_even[1+Nᵧ-n, n+1] = ExactReal(1)
                end
            end
            DF[:,:] = DγVtemp[:,:]*D_even[:,:]


        elseif j == 6
            DF[:,:] = project(v₂, codomain(DF))[:]
        elseif j ==7
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if n-2 >= 0
                                DF[:,:] = DF[:,:]  - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) ) ,domain(DF), codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j ==8
            DF[:,:] = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF))[:,:]) + project(Multiplication(λ*(v₂^0)),domain(DF),codomain(DF))[:,:] 
        elseif j ==9
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if m-2 >= 0
                                DF[:,:] = DF[:,:]  - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) ),domain(DF), codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 24
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if n-2 >= 0 && m-2>= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            elseif n-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ ,codomain(DF))[:]
                            elseif m-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ] * project(ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            end
                        end
                    end
                end
            end
        end
# i = 9 =============================================================================
    elseif i == 9
        if j == 6
            DF[:,:] = project(v₃, codomain(DF))[:]
        elseif j == 9
            DF[:,:]  = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF) )[:,:]) + project(Multiplication(λ*(v₃^0)), domain(DF),codomain(DF))[:,:] 

        elseif j == 10
            if eltype(α) ==  Interval{Float64}
                DF[:,:]  = -ω*interval.(Matrix(I, 2*order(codomain(DF))+1, 2*order(domain(DF))+1))
            else
                DF[:,:]  = -ω*Matrix(I, 2*order(codomain(DF))+1, 2*order(domain(DF))+1)
            end

        elseif j == 24
            DF[:,:] = project(-v₄,codomain(DF))[:]
        end
# i = 10 =============================================================================
    elseif i == 10
        if j == 1
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 && ℓ > 1
                            if n-2 >= 0 && m-2>= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₂[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            elseif n-2 >= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₂[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ ,codomain(DF))[:]
                            elseif m-2 >= 0
                                DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₂[n,m,ℓ] * project(ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            end
                        end
                    end
                end
            end
        elseif j == 2
                    
            codomainᵥ = codomain(DF)
            domainᵥ = domain(DF)
            Nᵥ =  order(codomainᵥ)
            Nᵧ =  order(domainᵥ)

            DγVtemp = zeros(eltype(DF), Fourier(order(domainᵥ), interval(1.0) ),codomainᵥ)


            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n-3 >= 0 
                                DγVtemp[:,:] = DγVtemp[:,:]  -ω*α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( ExactReal((n-1)*(n-2))*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:] 
                            end
                            if n-2 >= 0 && m-2 >= 0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:]
                            end

                        end
                    end
                end
            end

            D_even = zeros(eltype(DγVtemp) ,2*Nᵧ+1 , Nᵧ +1)
            for n = 0:Nᵧ
                if n == 0
                    D_even[1+Nᵧ, n+1] = ExactReal(1)
                else
                    D_even[1+Nᵧ+n, n+1] = ExactReal(1)
                    D_even[1+Nᵧ-n, n+1] = ExactReal(1)
                end
            end
            DF[:,:] = DγVtemp[:,:]*D_even[:,:]


        elseif j == 4
            codomainᵥ = codomain(DF)
            domainᵥ = domain(DF)
            Nᵥ =  order(codomainᵥ)
            Nᵧ =  order(domainᵥ)

            DγVtemp = zeros(eltype(DF), Fourier(order(domainᵥ), interval(1.0) ),codomainᵥ)

            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n-2 >= 0 && m-2 >= 0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:]  
                            end
                            if m-3 >=0
                                DγVtemp[:,:] = DγVtemp[:,:]  - ω*α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( ExactReal((m-1)*(m-2))*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ ) ,domain(DγVtemp), codomain(DγVtemp))[:,:] 
                            end
                        end
                    end
                end
            end
            D_even = zeros(eltype(DγVtemp) ,2*Nᵧ+1 , Nᵧ +1)
            for n = 0:Nᵧ
                if n == 0
                    D_even[1+Nᵧ, n+1] = ExactReal(1)
                else
                    D_even[1+Nᵧ+n, n+1] = ExactReal(1)
                    D_even[1+Nᵧ-n, n+1] = ExactReal(1)
                end
            end
            DF[:,:] = DγVtemp[:,:]*D_even[:,:]


        elseif j == 6
            DF[:,:] = project(v₄, codomain(DF))[:]
        elseif j == 7
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) ),domain(DF), codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 9
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if m-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) ),domain(DF), codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
        elseif j == 10
            DF[:,:] = Matrix(project(Derivative(1),domain(DF),codomain(DF), eltype(DF))[:,:]) + project(Multiplication(λ*(v₄^0)), domain(DF),codomain(DF))[:,:] 
        elseif j == 24
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n-2 >= 0 && m-2>= 0
                                DF[:,:]= DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃  ,codomain(DF))[:]
                            elseif n-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project(ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁  ,codomain(DF))[:]
                            elseif m-2 >= 0
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ] * project(ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ,codomain(DF))[:]
                            end
                        end
                    end
                end
            end
        end
# i = 11 =============================================================================
    elseif i == 11
        if j == 2
            Nw = order(codomain(DF))[1]
            Nt = order(domain(DF))[1]

            D_even = zeros(eltype(w),2*Nw+1 , Nt+1)
                for j = 0:min(Nw,Nt)
                    if j == 0
                        D_even[1+Nw, j+1] = ExactReal(1)
                    else
                        D_even[1+Nw+j, j+1] = ExactReal(1)
                        D_even[1+Nw-j, j+1] = ExactReal(1)
                    end
                end
                DF[(:,0),:] = -D_even

        elseif j == 6
            D₀₁ = Derivative(0,1)
            if eltype(α) ==  Interval{Float64}
                M = Sequence(Fourier(0,interval(1)) ⊗ Taylor(1), interval.([0, 1]))
            else
                M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])
            end
            DF[:,:] = project(M*(D₀₁*w₁),codomain(DF))[:]
            DF[(:,0) ,:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        elseif j == 7
            for k = -order(domain(DF)):order(domain(DF))
                for ℓ = -order(codomain(DF)[1]):order(codomain(DF)[1])
                    if ℓ == k
                        DF[(ℓ,1),k] = ExactReal(-1.0)
                    end
                end
            end
        elseif j == 11
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[:,:] = project(D₁₀,domain(DF),codomain(DF) )[:,:] 
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ*project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:]
           
            for i = -order(codomain(DF))[1]:order(codomain(DF))[1]
                for j = 0:1
                    for ℓ₁ = -order(domain(DF))[1]:order(domain(DF))[1]
                        for ℓ₂ = 0:order(domain(DF))[2]
                            if i == ℓ₁ && j == ℓ₂
                                DF[(i,j),(ℓ₁,ℓ₂)] = ExactReal(1)
                            else
                                DF[(i,j),(ℓ₁,ℓ₂)] = ExactReal(0)
                            end
                        end
                    end
                end
            end
        elseif j == 12
            if eltype(α) == Interval{Float64}
                DF[:,:] .= -ω*interval.(Matrix(I,size(DF[:,:] )) )
            else
                DF[:,:] .= -ω*Matrix(I,size(DF[:,:] )) 
            end
            DF[(:,0:1),(:,:)] .= ExactReal(0)
        elseif j == 24
            DF[:,:] = -project(w₂, codomain(DF))[:]
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        end
# i = 12 =============================================================================
    elseif i == 12
        if j == 1
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:] = DF[:,:]- ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₁[n,m,ℓ] * project(w₁^(n-1)*w₃^(m-1),codomain(DF))[:]
                        end
                    end
                end
            end
            DF[(:,0),1] .= ExactReal(0)
            DF[(:,1),1] .= ExactReal(0)
        elseif j == 3
            Nw = order(codomain(DF))[1]
            Nt = order(domain(DF))[1]

            D_odd= zeros(eltype(w),2*Nw +1,  Nt )
            if eltype(α) == Interval{Float64}
                for j = 1:min(Nw,Nt)
                    D_odd[Nw+j+1, j] = interval(-1im)
                    D_odd[Nw-j+1, j] = interval(1im)
                end
            else
                for j = 1:min(Nw,Nt)
                    D_odd[Nw+j+1, j] = -1im
                    D_odd[Nw-j+1, j] = 1im
                end
            end
             DF[(:,0),:] = -D_odd

        elseif j == 6
            D₀₁ = Derivative(0,1)
            if eltype(α) ==  Interval{Float64}
                M = Sequence(Fourier(0,interval(1)) ⊗ Taylor(1), interval.([0, 1]))
            else
                M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])
            end
            DF[:,:] = project(M*(D₀₁*w₂),codomain(DF))[:]
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        elseif j == 8
            for k = -order(domain(DF)):order(domain(DF))
                for ℓ = -order(codomain(DF)[1]):order(codomain(DF)[1])
                    if ℓ == k
                        DF[(ℓ,1),k] = ExactReal(-1.0)
                    end
                end
            end
        elseif j == 11 
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(w₁^(n-2)*w₃^(m-1)), domain(DF) ,codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(:,0:1),(:,:)] .= ExactReal(0)
        elseif j == 12
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[:,:] = project(D₁₀,domain(DF),codomain(DF) )[:,:] 
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ*project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:]
           
            for n = -order(codomain(DF))[1]:order(codomain(DF))[1]
                for m = 0:1
                    for ℓ₁ = -order(domain(DF))[1]:order(domain(DF))[1]
                        for ℓ₂ = 0:order(domain(DF))[2]
                            if n == ℓ₁ && m == ℓ₂
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal(1)
                            else
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal(0)
                            end
                        end
                    end
                end
            end
        elseif j == 13
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₁[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(w₁^(n-1)*w₃^(m-2)), domain(DF) ,codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(:,0:1),(:,:)] .= ExactReal(0)
        elseif j == 24
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ]*project(w₁^(n-1)*w₃^(m-1), codomain(DF))[:]
                        end
                    end
                end
            end
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        end
# i = 13 =============================================================================
    elseif i == 13
        if j == 4
            Nw = order(codomain(DF))[1]
            Nt = order(domain(DF))[1]

            D_even = zeros(eltype(w),2*Nw+1 , Nt+1)
                for m = 0:min(Nw,Nt)
                    if m == 0
                        D_even[1+Nw, m+1] = ExactReal(1)
                    else
                        D_even[1+Nw+m, m+1] = ExactReal(1)
                        D_even[1+Nw-m, m+1] = ExactReal(1)
                    end
                end
                DF[(:,0),:] = -D_even
        elseif j == 6
            D₀₁ = Derivative(0,1)
            if eltype(α) ==  Interval{Float64}
                M = Sequence(Fourier(0,interval(1)) ⊗ Taylor(1), interval.([0, 1]))
            else
                M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])
            end
            DF[:,:] = project(M*(D₀₁*w₃),codomain(DF))[:]
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        elseif j == 9
            for k = -order(domain(DF)):order(domain(DF))
                for ℓ = -order(codomain(DF)[1]):order(codomain(DF)[1])
                    if ℓ == k
                        DF[(ℓ,1),k] = ExactReal(-1.0)
                    end
                end
            end
        elseif j == 13
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[:,:] = project(D₁₀,domain(DF),codomain(DF) )[:,:] 
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ*project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:]
           
            for n = -order(codomain(DF))[1]:order(codomain(DF))[1]
                for m = 0:1
                    for ℓ₁ = -order(domain(DF))[1]:order(domain(DF))[1]
                        for ℓ₂ = 0:order(domain(DF))[2]
                            if n == ℓ₁ && m == ℓ₂
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal(1)
                            else
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal(0)
                            end
                        end
                    end
                end
            end
        elseif j == 14
            if eltype(α) == Interval{Float64}
                DF[:,:] .= -ω*interval.(Matrix(I,size(DF[:,:] )) )
            else
                DF[:,:] .= -ω*Matrix(I,size(DF[:,:] )) 
            end
            DF[(:,0:1),(:,:)] .= ExactReal(0)
        elseif j == 24
            DF[:,:] = -project(w₄, codomain(DF))[:]
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        end
# i = 14 =============================================================================
    elseif i == 14
        if j == 1
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*ω*M₂[n,m,ℓ] * project(w₁^(n-1)*w₃^(m-1),codomain(DF))[:]
                        end
                    end
                end
            end
            DF[(:,0),1] .= ExactReal( 0)
            DF[(:,1),1] .= ExactReal(0)
        elseif j == 5
            Nw = order(codomain(DF))[1]
            Nt = order(domain(DF))[1]

            D_odd= zeros(eltype(w),2*Nw +1,  Nt )
        
            if eltype(α) == Interval{Float64}
                for m = 1:min(Nw,Nt)
                    D_odd[Nw+m+1, m] = interval( -1im)
                     D_odd[Nw-m+1, m] = interval(1im)
               end
            else
                for m = 1:min(Nw,Nt)
                    D_odd[Nw+m+1, m] =-1im
                     D_odd[Nw-m+1, m] = 1im
               end
            end
             DF[(:,0),:] = -D_odd
             
        elseif j == 6
            D₀₁ = Derivative(0,1)
            if eltype(α) ==  Interval{Float64}
                M = Sequence(Fourier(0,interval(1)) ⊗ Taylor(1), interval.([0, 1]))
            else
                M = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])
            end
            DF[:,:] = project(M*(D₀₁*w₄),codomain(DF))[:]
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        elseif j == 10
            for k = -order(domain(DF)):order(domain(DF))
                for ℓ = -order(codomain(DF)[1]):order(codomain(DF)[1])
                    if ℓ == k
                        DF[(ℓ,1),k] = ExactReal(-1.0)
                    end
                end
            end
        elseif j == 11 
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(w₁^(n-2)*w₃^(m-1)), domain(DF) ,codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(:,0:1),(:,:)] .= ExactReal(0)
        elseif j == 13
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*ω*M₂[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(w₁^(n-1)*w₃^(m-2)), domain(DF) ,codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(:,0:1),(:,:)] .= ExactReal( 0)
        elseif j == 14
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[:,:] = project(D₁₀,domain(DF),codomain(DF) )[:,:] 
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ*project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:]
           
            for n = -order(codomain(DF))[1]:order(codomain(DF))[1]
                for m = 0:1
                    for ℓ₁ = -order(domain(DF))[1]:order(domain(DF))[1]
                        for ℓ₂ = 0:order(domain(DF))[2]
                            if n == ℓ₁ && m == ℓ₂
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal( 1)
                            else
                                DF[(n,m),(ℓ₁,ℓ₂)] = ExactReal(0)
                            end
                        end
                    end
                end
            end
        elseif j == 24
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ]*project(w₁^(n-1)*w₃^(m-1), codomain(DF))[:]
                        end
                    end
                end
            end
            DF[(:,0),:] .= ExactReal(0)
            DF[(:,1),:] .= ExactReal(0)
        end
# i = 15 =============================================================================
    elseif i == 15
        if j == 18
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = (D₁₀*p₁)( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) + (D₀₁*p₁)( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )
            else
                DF[1,1] = (D₁₀*p₁)( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*p₁)( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
            end

        elseif j == 19
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = interval(1im)*(D₁₀*p₁)( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) - interval(1im)*(D₀₁*p₁)( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )    
            else
                DF[1,1] = 1im*(D₁₀*p₁)( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*p₁)( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
            end
        elseif j == 20
            DF[:,:] .= ExactReal(-2)
            DF[:,0] .= ExactReal(-1)
        elseif j == 26
            Nₚ = order(domain(DF))
            n = ExactReal.(reshape(repeat(0:Nₚ[1],1,Nₚ[2]+1), (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            m = ExactReal.(reshape(repeat(transpose(0:Nₚ[2]),Nₚ[1]+1,1) , (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            if eltype(α) == Interval{Float64}
                DF[:,:] = transpose((σ₁[1] + interval(1im)*σ₂[1] ).^n.*(σ₁[1] - interval(1im)*σ₂[1]).^m)
            else
                DF[:,:] = transpose((σ₁[1] + 1im*σ₂[1] ).^n.*(σ₁[1] - 1im*σ₂[1]).^m)
            end
        end
# i = 16 =============================================================================
    elseif i == 16
        if j == 18
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = (D₁₀*component(p,2))( σ₁[1] +  interval(1im)*σ₂[1], σ₁[1] -  interval(1im)*σ₂[1] ) + (D₀₁*component(p,2))( σ₁[1] +  interval(1im)*σ₂[1], σ₁[1] -  interval(1im)*σ₂[1] )
            else
                DF[1,1] = (D₁₀*component(p,2))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*component(p,2))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
            end

        elseif j == 19
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] =  interval(1im)*(D₁₀*component(p,2))( σ₁[1] +  interval(1im)*σ₂[1], σ₁[1] -  interval(1im)*σ₂[1] ) -  interval(1im)*(D₀₁*component(p,2))( σ₁[1] +  interval(1im)*σ₂[1], σ₁[1] -  interval(1im)*σ₂[1] )    
            else
                DF[1,1] = 1im*(D₁₀*component(p,2))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*component(p,2))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
            end

        elseif j == 21
            DF[:,:] .= ExactReal(-2)
            DF[:,0] .= ExactReal(-1)
        elseif j == 27
            Nₚ = order(domain(DF))
            n = ExactReal.(reshape(repeat(0:Nₚ[1],1,Nₚ[2]+1), (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            m = ExactReal.(reshape(repeat(transpose(0:Nₚ[2]),Nₚ[1]+1,1) , (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            if eltype(α) == Interval{Float64}
                DF[:,:] = transpose((σ₁[1] +  interval(1im)*σ₂[1] ).^n.*(σ₁[1] -  interval(1im)*σ₂[1]).^m)
            else
                DF[:,:] = transpose((σ₁[1] + 1im*σ₂[1] ).^n.*(σ₁[1] - 1im*σ₂[1]).^m)
            end
        end
# i = 17 =============================================================================
    elseif i == 17
        if j == 18
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = (D₁₀*component(p,3))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) + (D₀₁*component(p,3))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )
            else
                DF[1,1] = (D₁₀*component(p,3))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*component(p,3))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
            end
        elseif j == 19
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = interval(1im)*(D₁₀*component(p,3))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) - interval(1im)*(D₀₁*component(p,3))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )    
            else
                DF[1,1] = 1im*(D₁₀*component(p,3))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*component(p,3))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
            end

        elseif j == 22
            DF[:,:] .= ExactReal(-2)
            DF[:,0] .= ExactReal(-1)
        elseif j == 28
            Nₚ = order(domain(DF))
            n = ExactReal.(reshape(repeat(0:Nₚ[1],1,Nₚ[2]+1), (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            m = ExactReal.(reshape(repeat(transpose(0:Nₚ[2]),Nₚ[1]+1,1) , (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            if eltype(α) == Interval{Float64}
                DF[:,:] = transpose((σ₁[1] + interval(1im)*σ₂[1] ).^n.*(σ₁[1] - interval(1im)*σ₂[1]).^m)
            else
                DF[:,:] = transpose((σ₁[1] + 1im*σ₂[1] ).^n.*(σ₁[1] - 1im*σ₂[1]).^m)
            end
        end
# i = 18 =============================================================================
    elseif i == 18
        if j == 18
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = (D₁₀*component(p,4))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) + (D₀₁*component(p,4))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )
            else
                DF[1,1] = (D₁₀*component(p,4))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*component(p,4))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
            end

        elseif j == 19
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                DF[1,1] = interval(1im)*(D₁₀*component(p,4))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] ) - interval(1im)*(D₀₁*component(p,4))( σ₁[1] + interval(1im)*σ₂[1], σ₁[1] - interval(1im)*σ₂[1] )    
            else
                DF[1,1] = 1im*(D₁₀*component(p,4))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*component(p,4))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
            end

        elseif j == 23
            DF[:,:] .= ExactReal(-2)
            DF[:,0] .= ExactReal(-1)
        elseif j == 29
            Nₚ = order(domain(DF))
            n = ExactReal.(reshape(repeat(0:Nₚ[1],1,Nₚ[2]+1), (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            m = ExactReal.(reshape(repeat(transpose(0:Nₚ[2]),Nₚ[1]+1,1) , (Nₚ[1]+1)*(Nₚ[2]+1)  ,1)[:])
            if eltype(α) == Interval{Float64}
                DF[:,:] = transpose((σ₁[1] + interval(1im)*σ₂[1] ).^n.*(σ₁[1] - interval(1im)*σ₂[1]).^m)
            else
                DF[:,:] = transpose((σ₁[1] + 1im*σ₂[1] ).^n.*(σ₁[1] - 1im*σ₂[1]).^m)
            end
        end
# i = 19 =============================================================================
    elseif i == 19
        if j == 11
            N_w = order(domain(DF))   
            wᵢ = project(component(w,1) , domain(DF))
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = (interval(-0.95).^m[(0,:)])[:]
            else
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)])[:]
                 DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)])[:]
                 DF[0,(0,:)] = ((-0.95).^m[(0,:)])[:]
            end
        elseif j == 15
            Ψ = ExactReal(1/2)*project( a₂ ,Chebyshev(order(codomain(DF))+1) )
            DF[1:end,:] = DF[1:end,:] + Ψ[2:end]  - Ψ[0:end-2] 

        elseif j == 16
            N_w = order(w)[1] 
            wᵢ = component(w,1)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)] .- interval(1))).*interval(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
            
        elseif j == 17
            N_w = order(w)[1] 
            wᵢ = component(w,1)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(interval(1im)*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)]  - interval(1im)*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
                        
        elseif j == 20
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            for n = 0:Nₐ_domain
                 for m = 0:Nₐ_codomain
                    if n == m
                        DF[ m, n] = ExactReal(2*n)
                    end
                end
            end            
            IC = 2*ones(1,Nₐ_domain+1) 
            IC[1] = 1
            IC = -IC.*((-1).^transpose(0:Nₐ_domain))
            DF[0,0:Nₐ_domain ] = ExactReal.(IC)
        elseif j == 21
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(DaΨ[:,:],1)
                for m in axes(DaΨ[:,:],2)
                    if n == m
                        DaΨ[n-1,m-1] = L/ExactReal(2)
                    end
                end
            end
            T = zeros(eltype(w), domain(DaΨ), domain(DaΨ) )
            for k = 0:order(domain(DaΨ))
                for m = 1:order(domain(DaΨ))
                    if k == m-1 
                        T[m,k] = ExactReal(-1.)
                    end
                    if k == m+1 
                        T[m,k] = ExactReal(1.)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]
        

        end
# i = 20 =============================================================================
    elseif i == 20
        if j == 1
            DαΨ = Sequence(Chebyshev(order( codomain(DF) )+1),zeros(eltype(x),order(codomain(DF) )+2))
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 && ℓ > 1
                            DαΨ = DαΨ + (L/ExactReal(2))*(ℓ-1)*α^ExactReal(ℓ-2)*M₁[n,m,ℓ] * project(a₁^(n-1)*a₃^(m-1),space(DαΨ))
                        end
                    end
                end
            end
            DF[1:end,1] = DαΨ[2:end] - DαΨ[0:end-2] 
            DF[0,1] = ExactReal(0)
        elseif j == 12
            N_w = order(domain(DF))   
            wᵢ = project(component(w,2) , domain(DF))
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = (interval(-0.95).^m[(0,:)])[:]
            else
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = ((-0.95).^m[(0,:)])[:]
            end
        elseif j == 15
            Ψ = Sequence(Chebyshev(order(codomain(DF))+1),zeros(eltype(α), order(codomain(DF))+2))
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ  in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            Ψ = Ψ + ExactReal(1/2)*M₁[n,m,ℓ] * α^ExactReal(ℓ-1)*project(a₁^(n-1)*a₃^(m-1),space(Ψ))
                        end
                    end
                end
            end
            DF[1:end,:] = DF[1:end,:] + Ψ[2:end]  - Ψ[0:end-2] 
        elseif j == 16
            N_w = order(w)[2] 
            wᵢ = component(w,2)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
            
        elseif j == 17
            N_w = order(w)[2] 
            wᵢ = component(w,2)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(interval(1im)*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)]  - interval(1im)*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .- interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
        elseif j == 20
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if n >= 2
                            DaΨ[:,:] = DaΨ[:,:] + α^ExactReal(ℓ-1)*L/ExactReal(2)*M₁[n,m,ℓ] * project(  Multiplication(ExactReal(n-1)*a₁^(n-2)*a₃^(m-1))   ,Chebyshev(Nₐ_domain+1),Chebyshev(Nₐ_codomain +1 ))[:,:]
                            end
                        end
                    end
                end
            end
            T = (zeros(eltype(α), domain(DaΨ), domain(DaΨ) ))
            for k = 0:order(domain(DaΨ))
                for m = 1:order(domain(DaΨ))
                    if k == m-1 
                        T[m,k] = ExactReal(-1.0)
                    end
                    if k == m+1 
                        T[m,k] = ExactReal(1.0)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]

        elseif j == 21
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            for n = 0:Nₐ_domain
                 for m = 0:Nₐ_codomain
                    if n == m
                        DF[ m, n] = ExactReal(2*n)
                    end
                end
            end            
            IC = 2*ones(1,Nₐ_domain+1) 
            IC[1] = 1
            IC = -IC.*((-1).^transpose(0:Nₐ_domain))
            DF[0,0:Nₐ_domain ] = ExactReal.(IC)
        elseif j == 22
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0
                            if m>= 2
                                DaΨ[:,:] = DaΨ[:,:] + α^ExactReal(ℓ-1)*L/ExactReal(2)*M₁[n,m,ℓ] * project(  Multiplication(ExactReal(m-1)*a₁^(n-1)*a₃^(m-2))   ,Chebyshev(Nₐ_domain+1),Chebyshev(Nₐ_codomain +1 ))[:,:]
                            end
                        end
                    end
                end
            end
            T = (zeros( eltype(α) ,domain(DaΨ), domain(DaΨ) ))
            for k = 0:order(domain(DaΨ))
                for m = 1:order(domain(DaΨ))
                    if k == m-1 
                        T[m,k] = ExactReal(-1.)
                    end
                    if k == m+1 
                        T[m,k] = ExactReal(1.)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]
        end
# i = 21 =============================================================================
    elseif i == 21
        if j == 13
            N_w = order(domain(DF))   
            wᵢ = project(component(w,3) , domain(DF))
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = (interval(-0.95).^m[(0,:)])[:]
            else
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = ((-0.95).^m[(0,:)])[:]
            end
        elseif j == 15
            Ψ = ExactReal(1/2)*project( a₄ ,Chebyshev(order(codomain(DF))+1) )
            DF[1:end,:] = DF[1:end,:] + Ψ[2:end]  - Ψ[0:end-2] 
        elseif j == 16
            N_w = order(w)[3] 
            wᵢ = component(w,3)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
            
        elseif j == 17
            N_w = order(w)[3] 
            wᵢ = component(w,3)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(interval(1im)*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)]  - interval(1im)*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
        elseif j == 22
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            for n = 0:Nₐ_domain
                 for m = 0:Nₐ_codomain
                    if n == m
                        DF[ m, n] = ExactReal(2*n)
                    end
                end
            end            
            IC = 2*ones(1,Nₐ_domain+1) 
            IC[1] = 1
            IC = -IC.*((-1).^transpose(0:Nₐ_domain))
            DF[0,0:Nₐ_domain ] =  ExactReal.(IC)
        elseif j == 23
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(DaΨ[:,:],1)
                for m in axes(DaΨ[:,:],2)
                    if n == m
                        DaΨ[n-1,m-1] = L/ExactReal(2)
                    end
                end
            end
            T = zeros( eltype(α), codomain(DaΨ), codomain(DaΨ) )
            for k = 0:order(codomain(DaΨ))
                for m = 1:order(codomain(DaΨ))
                    if k == m-1 
                        T[m,k] =  ExactReal(-1)
                    end
                    if k == m+1 
                        T[m,k] =  ExactReal(1)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]
        end
# i = 22 =============================================================================
    elseif i == 22
        if j == 1
            DαΨ = Sequence(Chebyshev(order( codomain(DF) )+1),zeros(eltype(x),order(codomain(DF) )+2))
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 && ℓ > 1
                            DαΨ = DαΨ + (L/ExactReal(2))*ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*M₂[n,m,ℓ] * project(a₁^(n-1)*a₃^(m-1),space(DαΨ))
                        end
                    end
                end
            end
            DF[1:end,1] = DαΨ[2:end] - DαΨ[0:end-2] 
            DF[0,1] = ExactReal(0)
        elseif j == 14
            N_w = order(domain(DF))   
            wᵢ = project(component(w,4) , domain(DF))
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + interval(1im)*θ₂[1]).^n[(1:N_w[1],:)]).*interval(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - interval(1im)*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*interval(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = (interval(-0.95).^m[(0,:)])[:]
            else
                DF[0,((1:N_w[1]),:)] = transpose(((θ₁[1] + 1im*θ₂[1]).^n[(1:N_w[1],:)]).*(-0.95).^m[(1:N_w[1],:)])[:]
                DF[0,(-N_w[1]:-1,:)] = transpose(((θ₁[1] - 1im*θ₂[1]).^abs.(n[(-N_w[1]:-1,:)])).*(-0.95).^m[(-N_w[1]:-1,:)])[:]
                DF[0,(0,:)] = ((-0.95).^m[(0,:)])[:]
            end
        elseif j == 15
            Ψ = Sequence(Chebyshev(order(codomain(DF))+1),zeros(eltype(α),order(codomain(DF))+2))
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ  in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            Ψ = Ψ + ExactReal(1/2)*M₂[n,m,ℓ] * α^ExactReal(ℓ-1)*project(a₁^(n-1)*a₃^(m-1),space(Ψ))
                        end
                    end
                end
            end
            DF[1:end,:] = DF[1:end,:] + Ψ[2:end]  - Ψ[0:end-2] 
        elseif j == 16
            N_w = order(w)[4] 
            wᵢ = component(w,4)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
                    
        elseif j == 17
            N_w = order(w)[4] 
            wᵢ = component(w,4)
            n = ExactReal.(Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            m = ExactReal.(Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:]))
            if eltype(α) == Interval{Float64}
                DF[0,1] = sum(interval(1im)*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + interval(1im)*θ₂[1]).^(n[(1:N_w[1],:)].-interval(1))).*interval(-0.95).^m[(1:N_w[1],:)]  - interval(1im)*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - interval(1im)*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-interval(1))).*interval(-0.95).^m[(-N_w[1]:-1,:)] ) 
            else
                DF[0,1] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
            end
        elseif j == 20
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if n >= 2
                                DaΨ[:,:] = DaΨ[:,:] + α^ExactReal(ℓ-1)*L/ExactReal(2)*M₂[n,m,ℓ] * project(  Multiplication(ExactReal(n-1)*a₁^(n-2)*a₃^(m-1))   ,Chebyshev(Nₐ_domain+1),Chebyshev(Nₐ_codomain +1 ))[:,:]
                            end
                        end
                    end
                end
            end
            T = zeros( eltype(α),  domain(DaΨ), domain(DaΨ) )
            for k = 0:order(domain(DaΨ))
                for m = 1:order(domain(DaΨ))
                    if k == m-1 
                        T[m,k] = ExactReal(-1)
                    end
                    if k == m+1 
                        T[m,k] = ExactReal(1)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]
        elseif j == 22
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            DaΨ = zeros( eltype(w), Chebyshev(Nₐ_domain+1), Chebyshev(Nₐ_codomain +1 )  )
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0
                            if m>= 2
                                DaΨ[:,:] = DaΨ[:,:] + α^ExactReal(ℓ-1)*L/ExactReal(2)*M₂[n,m,ℓ] * project(  Multiplication(ExactReal(m-1)*a₁^(n-1)*a₃^(m-2))   ,Chebyshev(Nₐ_domain+1),Chebyshev(Nₐ_codomain +1 ))[:,:]
                            end
                        end
                    end
                end
            end
            T = zeros(eltype(α),domain(DaΨ), domain(DaΨ) )
            for k = 0:order(domain(DaΨ))
                for m = 1:order(domain(DaΨ))
                    if k == m-1 
                        T[m,k] = ExactReal(-1)
                    end
                    if k == m+1 
                        T[m,k] = ExactReal(1)
                    end
                end
            end
            DF[:,:] = project(T*DaΨ, domain(DF) ,codomain(DF))[:,:]
        elseif j == 23
            Nₐ_domain = order(domain(DF))
            Nₐ_codomain = order(codomain(DF))
            for n = 0:Nₐ_domain
                 for m = 0:Nₐ_codomain
                    if n == m
                        DF[ m, n] = ExactReal(2*n)
                    end
                end
            end            
            IC = 2*ones(1,Nₐ_domain+1) 
            IC[1] = 1
            IC = -IC.*((-1).^transpose(0:Nₐ_domain))
            DF[0,0:Nₐ_domain ] = ExactReal.(IC)
        end
# i = 23 =============================================================================
    elseif i == 23
        if j == 18
            DF[1,1] .= ExactReal(-2)*σ₁
        elseif j == 19
            DF[1,1] .= ExactReal(-2)*σ₂
        end
# i = 24 =============================================================================
    elseif i == 24
        if j == 16
            DF[1,1] = ExactReal(-2)*θ₁[1]
        elseif j == 17
            DF[1,1] = ExactReal(-2)*θ₂[1]
        end
# i = 25 =============================================================================
    elseif i == 25
        if j == 1
                    
            M = zeros(eltype(x),4,4)
            for ℓ in axes(M₁,3)
                if ℓ > 1 
                    if size(M₁,1) > 1
                        M[2,1] = M[2,1] + M₁[2,1,ℓ]*ExactReal(ℓ-1)*α^ExactReal(ℓ-2)
                    end
                    if size(M₁,2) > 1
                        M[2,3] = M[2,3] + M₁[1,2,ℓ]*ExactReal(ℓ-1)*α^ExactReal(ℓ-2)
                    end
                end
            end
            for ℓ in axes(M₂,3)
                if ℓ > 1 
                    if size(M₂,1) > 1
                        M[4,1] = M[4,1] + M₂[2,1,ℓ]*ExactReal(ℓ-1)*α^ExactReal(ℓ-2)
                    end
                    if size(M₂,2) > 1
                        M[4,3] = M[4,3] + M₂[1,2,ℓ]*ExactReal(ℓ-1)*α^ExactReal(ℓ-2)
                    end
                end
            end
            DF[:,:] =  [ ExactReal(0);- M*component(eigenpairs,2:5)[:]; ExactReal(0); - M*component(eigenpairs,7:10)[:] ]
        elseif j == 25
            DF = DE!(DF, eigenpairs, M₁, M₂,α, star)
        end
# i = 26 =============================================================================
    elseif i == 26
        if j == 25
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), interval.([0, 1]))
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), interval.([0, 1]))
            else
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])
            end

            Dλ₁P = project(M₁₀*(D₁₀*p₁), codomain(DF) )
            Dλ₁P[(0,0)] = ExactReal(0)
            Dλ₁P[(1,0)] = ExactReal(0)
            Dλ₁P[(0,1)] = ExactReal(0)
            component(DF,1)[:,:] = Dλ₁P[:]

            Dλ₂P = project(M₀₁*(D₀₁*p₁), codomain(DF) )
            Dλ₂P[(0,0)] = ExactReal(0)
            Dλ₂P[(1,0)] = ExactReal(0)
            Dλ₂P[(0,1)] = ExactReal(0)
            component(DF,6)[:,:] = Dλ₂P[:]

            component(DF,2)[(1,0),:] .= ExactReal(-1)
            component(DF,7)[(0,1),:] .= ExactReal(-1) 
    
        elseif j == 26
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ₂*ExactReal.(project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:])
            DF[(1:order( DF )[2][1],:),(:,:)] = DF[(1:order( DF )[2][1],:),(:,:)] + λ₁*ExactReal.(project(D₁₀,domain(DF),Taylor( order(DF)[2][1]-1  ) ⊗ codomain(DF)[2]  )[:,:])
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(0,0),(0,0)] = ExactReal(1)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(1,0)] = ExactReal(1)
            DF[(0,1),(:,:)] .= ExactReal(0)
            DF[(0,1),(0,1)] = ExactReal(1)
        elseif j == 27
            if eltype(α) == Interval{Float64}
                DF[:,:] .= interval.( -1.0*Matrix(I,size(DF[:,:])))
            else
                DF[:,:] .= -1.0*Matrix(I,size(DF[:,:]))
            end
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        end
# i = 27 =============================================================================
    elseif i == 27
        if j == 1
            for n  in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*M₁[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),codomain(DF))[:]

                        end
                    end
                end
            end
            DF[(0,0),1] = ExactReal(0)
            DF[(1,0),1] = ExactReal(0)
            DF[(0,1),1] = ExactReal(0)
        
        elseif j == 25
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)

            if eltype(α) == Interval{Float64}
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), interval.([0, 1]))
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), interval.([0, 1]))
            else
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])
            end

            Dλ₁P = project(M₁₀*(D₁₀*p₂), codomain(DF) )
            Dλ₁P[(0,0)] = ExactReal(0)
            Dλ₁P[(1,0)] = ExactReal(0)
            Dλ₁P[(0,1)] = ExactReal(0)
            component(DF,1)[:,:] = Dλ₁P[:]

            Dλ₂P = project(M₀₁*(D₀₁*p₂), codomain(DF) )
            Dλ₂P[(0,0)] = ExactReal(0)
            Dλ₂P[(1,0)] = ExactReal(0)
            Dλ₂P[(0,1)] = ExactReal(0)
            component(DF,6)[:,:] = Dλ₂P[:]

            component(DF,3)[(1,0),:] .=ExactReal( -1)
            component(DF,8)[(0,1),:] .= ExactReal(-1 )
        elseif j == 26
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(p₁^(n-2)*p₃^(m-1)), domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        elseif j == 27
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ₂*ExactReal.(project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:])
            DF[(1:order( DF )[2][1],:),(:,:)] = DF[(1:order( DF )[2][1],:),(:,:)] + λ₁*ExactReal.(project(D₁₀,domain(DF),Taylor( order(DF)[2][1]-1  ) ⊗ codomain(DF)[2]  )[:,:])
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(0,0),(0,0)] = ExactReal(1)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(1,0)] = ExactReal(1)
            DF[(0,1),(:,:)] .= ExactReal(0)
            DF[(0,1),(0,1)] = ExactReal(1)
        elseif j == 28
            for n in axes(M₁,1)
                for m in axes(M₁,2)
                    for ℓ in axes(M₁,3)
                        if M₁[n,m,ℓ] != 0 
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₁[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(p₁^(n-1)*p₃^(m-2)), domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        end
# i = 28 =============================================================================
    elseif i == 28
        if j == 25
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            if eltype(α) == Interval{Float64}
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), interval.([0, 1]))
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), interval.([0, 1]))
            else
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])
            end

            Dλ₁P = project(M₁₀*(D₁₀*p₃), codomain(DF) )
            Dλ₁P[(0,0)] = ExactReal(0)
            Dλ₁P[(1,0)] = ExactReal(0)
            Dλ₁P[(0,1)] = ExactReal(0)
            component(DF,1)[:,:] = Dλ₁P[:]

            Dλ₂P = project(M₀₁*(D₀₁*p₃), codomain(DF) )
            Dλ₂P[(0,0)] = ExactReal(0)
            Dλ₂P[(1,0)] = ExactReal(0)
            Dλ₂P[(0,1)] = ExactReal(0)
            component(DF,6)[:,:] = Dλ₂P[:]

            component(DF,4)[(1,0),:] .= ExactReal(-1)
            component(DF,9)[(0,1),:] .= ExactReal(-1) 
        elseif j == 28
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ₂*ExactReal.(project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:])
            DF[(1:order( DF )[2][1],:),(:,:)] = DF[(1:order( DF )[2][1],:),(:,:)] + λ₁*ExactReal.(project(D₁₀,domain(DF),Taylor( order(DF)[2][1]-1  ) ⊗ codomain(DF)[2]  )[:,:])
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(0,0),(0,0)] = ExactReal(1)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(1,0)] = ExactReal(1)
            DF[(0,1),(:,:)] .= ExactReal(0)
            DF[(0,1),(0,1)] = ExactReal(1)
        elseif j == 29
            if eltype(α) == Interval{Float64}
                DF[:,:] .= interval.( -1.0*Matrix(I,size(DF[:,:])))
            else
                DF[:,:] .= -1.0*Matrix(I,size(DF[:,:]))
            end
            DF[(0,0),(:,:)] .=ExactReal( 0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        end
# i = 29 =============================================================================
    elseif i == 29
        if j == 1
            for n  in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 && ℓ > 1
                            DF[:,:] = DF[:,:] - ExactReal(ℓ-1)*α^ExactReal(ℓ-2)*M₂[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),codomain(DF))[:]
                        end
                    end
                end
            end
            DF[(0,0),1] = ExactReal(0)
            DF[(1,0),1] = ExactReal(0)
            DF[(0,1),1] = ExactReal(0)
        elseif j == 25
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)

            if eltype(α) == Interval{Float64}
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), interval.([0, 1]))
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), interval.([0, 1]))
            else
                M₁₀ = Sequence(Taylor(1) ⊗ Taylor(0), [0, 1])
                M₀₁ = Sequence(Taylor(0) ⊗ Taylor(1), [0, 1])
            end

            Dλ₁P = project(M₁₀*(D₁₀*p₄), codomain(DF) )
            Dλ₁P[(0,0)] = ExactReal(0)
            Dλ₁P[(1,0)] = ExactReal(0)
            Dλ₁P[(0,1)] = ExactReal(0)
            component(DF,1)[:,:] = Dλ₁P[:]

            Dλ₂P = project(M₀₁*(D₀₁*p₄), codomain(DF) )
            Dλ₂P[(0,0)] = ExactReal(0)
            Dλ₂P[(1,0)] = ExactReal(0)
            Dλ₂P[(0,1)] = ExactReal(0)
            component(DF,6)[:,:] = Dλ₂P[:]

            component(DF,5)[(1,0),:] .= ExactReal(-1)
            component(DF,10)[(0,1),:] .= ExactReal(-1) 
        elseif j == 26
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 
                            if n >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ]*ExactReal(n-1)*project( Multiplication(p₁^(n-2)*p₃^(m-1)), domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        elseif j == 28
            for n in axes(M₂,1)
                for m in axes(M₂,2)
                    for ℓ in axes(M₂,3)
                        if M₂[n,m,ℓ] != 0 
                            if m >= 2
                                DF[:,:] = DF[:,:] - α^ExactReal(ℓ-1)*M₂[n,m,ℓ]*ExactReal(m-1)*project( Multiplication(p₁^(n-1)*p₃^(m-2)), domain(DF),codomain(DF))[:,:]
                            end
                        end
                    end
                end
            end
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(0,1),(:,:)] .= ExactReal(0)
        elseif j == 29
            D₁₀ = Derivative(1,0)
            D₀₁ = Derivative(0,1)
            DF[(:,1:order( DF )[2][2]),(:,:)] = DF[(:,1:order( DF )[2][2]),(:,:)] + λ₂*ExactReal.(project(D₀₁,domain(DF),codomain(DF)[1] ⊗ Taylor( order(DF)[2][2]-1  )  )[:,:])
            DF[(1:order( DF )[2][1],:),(:,:)] = DF[(1:order( DF )[2][1],:),(:,:)] + λ₁*ExactReal.(project(D₁₀,domain(DF),Taylor( order(DF)[2][1]-1  ) ⊗ codomain(DF)[2]  )[:,:])
            DF[(0,0),(:,:)] .= ExactReal(0)
            DF[(0,0),(0,0)] = ExactReal(1)
            DF[(1,0),(:,:)] .= ExactReal(0)
            DF[(1,0),(1,0)] = ExactReal(1)
            DF[(0,1),(:,:)] .= ExactReal(0)
            DF[(0,1),(0,1)] = ExactReal(1)
        end
    end

    return DF
end

function DEP(DF, α, eigenpairs, p , star, M₁, M₂)
    DE = component(DF , 1 , 1)
    DeigenpairsP = component(DF , 2:5 , 1)
    DₚP = component(DF , 2:5 , 2:5)
    DE = DE!(DE, eigenpairs, M₁, M₂,α, star)
    DeigenpairsP, DₚP = DP!(DeigenpairsP, DₚP,p, eigenpairs , M₁, M₂, α)
    return DF
end

function DαEP(DαF, α, eigenpairs, p , M₁, M₂)
    p₁, p₂, p₃, p₄ = eachcomponent(p)
    

    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0 && ℓ > 1
                    component(DαF,3)[:] = component(DαF,3)[:] - (ℓ-1)*α^(ℓ-2)*M₁[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),space(component(DαF,3)))[:]
                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0 && ℓ > 1
                    component(DαF,5)[:] = component(DαF,5)[:] - (ℓ-1)*α^(ℓ-2)*M₂[n,m,ℓ] * project(p₁^(n-1)*p₃^(m-1),space(component(DαF,5)))[:]
                end
            end
        end
    end
    component(DαF,3)[(0,0)] = 0
    component(DαF,3)[(1,0)] = 0
    component(DαF,3)[(0,1)] = 0
    component(DαF,5)[(0,0)] = 0
    component(DαF,5)[(1,0)] = 0
    component(DαF,5)[(0,1)] = 0

    
    M = zeros(eltype(eigenpairs),4,4)
    for ℓ in axes(M₁,3)
        if ℓ > 1 
            if size(M₁,1) > 1
                M[2,1] = M[2,1] + M₁[2,1,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
            if size(M₁,2) > 1
                M[2,3] = M[2,3] + M₁[1,2,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
        end
    end
    for ℓ in axes(M₂,3)
        if ℓ > 1 
            if size(M₂,1) > 1
                M[4,1] = M[4,1] + M₂[2,1,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
            if size(M₂,2) > 1
                M[4,3] = M[4,3] + M₂[1,2,ℓ]*(ℓ-1)*α^(ℓ-2)
            end
        end
    end
    component(DαF,1)[:] =  [ 0;- M*component(eigenpairs,2:5)[:]; 0; - M*component(eigenpairs,7:10)[:] ]

    return DαF
end

function DV2!(  DλV, DᵥV,DωV,DᵧV, DλC, DᵥC,DωC,DᵧC,Dc₀C , DcC,  ω, x, γ, M₁, M₂,α,  κ)
    # Parameters 
    λ = component(x,1)[1]
    v = component(x,2:5)
    c₀ = component(x,6)[1]
    c = component(x,7:10)
        γ_ext =  γ_extened_orientable!(γ)
        v₁,v₂,v₃,v₄ = eachcomponent(v)
        c₁,c₂,c₃,c₄ = eachcomponent(c)

        γ₁ = component(γ_ext,1)
        γ₃ = component(γ_ext,3)
        D = Derivative(1)
        
    # DωV 
        component(DωV,1)[1] = 0
        component(DωV,2)[:] = project(-v₂,space(DωV)[2])[:]
        component(DωV,3)[:] .= 0
        component(DωV,4)[:] = project(-v₄,space(DωV)[4])[:]
        component(DωV,5)[:] .= 0
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-2 >= 0 && m-2>= 0
                            component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[3])[:]
                        elseif n-2 >= 0
                            component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ ,space(DωV)[3])[:]
                        elseif m-2 >= 0
                            component(DωV,3)[:] = component(DωV,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[3])[:]
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
                            component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃  ,space(DωV)[5])[:]
                        elseif n-2 >= 0
                            component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*v₁  ,space(DωV)[5])[:]
                        elseif m-2 >= 0
                            component(DωV,5)[:] = component(DωV,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*v₃ ,space(DωV)[5])[:]
                        end
                    end
                end
            end
        end
    
    # DλV
        component(DλV,1)[1]= 0.0
        component(DλV,2)[:] = project(v₁, space(DλV)[2] )[:]
        component(DλV,3)[:] = project(v₂, space(DλV)[3] )[:]
        component(DλV,4)[:] = project(v₃, space(DλV)[4] )[:]
        component(DλV,5)[:] = project(v₄, space(DλV)[5] )[:]
     
    # DᵥV
        component(DᵥV,2,1)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[1],codomain(DᵥV)[2], eltype(DᵥV) )[:,:]) + project(Multiplication(λ*(v₁^0)),domain(DᵥV)[1],codomain(DᵥV)[2])[:,:] 
        component(DᵥV,2,2)[:,:] = -ω*Matrix(I, 2*order(codomain(DᵥV)[2])+1, 2*order(domain(DᵥV)[2])+1)
        component(DᵥV,2,3)[:,:] = zeros(2*order(codomain(DᵥV)[2])+1,2*order(domain(DᵥV)[3])+1)
        component(DᵥV,2,4)[:,:] = zeros(2*order(codomain(DᵥV)[2])+1,2*order(domain(DᵥV)[4])+1)
    
        component(DᵥV,3,1)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[1])+1)
        component(DᵥV,3,2)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[2],codomain(DᵥV)[3], eltype(DᵥV))[:,:]) + project(Multiplication(λ*(v₂^0)),domain(DᵥV)[2],codomain(DᵥV)[3])[:,:] 
        component(DᵥV,3,3)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[3])+1)
        component(DᵥV,3,4)[:,:] = zeros(2*order(codomain(DᵥV)[3])+1,2*order(domain(DᵥV)[4])+1)
    
        component(DᵥV,4,1)[:,:] = zeros(2*order(codomain(DᵥV)[4])+1,2*order(domain(DᵥV)[1])+1)
        component(DᵥV,4,2)[:,:] = zeros(2*order(codomain(DᵥV)[4])+1,2*order(domain(DᵥV)[2])+1)
        component(DᵥV,4,3)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[3],codomain(DᵥV)[4], eltype(DᵥV) )[:,:]) + project(Multiplication(λ*(v₃^0)), domain(DᵥV)[3],codomain(DᵥV)[4])[:,:] 
        component(DᵥV,4,4)[:,:] = -ω*Matrix(I, 2*order(codomain(DᵥV)[4])+1, 2*order(domain(DᵥV)[4])+1)
    
        component(DᵥV,5,1)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[1])+1)
        component(DᵥV,5,2)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[2])+1)
        component(DᵥV,5,3)[:,:] = zeros(2*order(codomain(DᵥV)[5])+1,2*order(domain(DᵥV)[3])+1)
        component(DᵥV,5,4)[:,:] = Matrix(project(Derivative(1),domain(DᵥV)[4],codomain(DᵥV)[5], eltype(DᵥV))[:,:]) + project(Multiplication(λ*(v₄^0)), domain(DᵥV)[4],codomain(DᵥV)[5])[:,:] 
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-2 >= 0
                            component(DᵥV,3,1)[:,:] = component(DᵥV,3,1)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ) ,domain(DᵥV)[1], codomain(DᵥV)[3])[:,:]
                        end
                        if m-2 >= 0
                            component(DᵥV,3,3)[:,:] = component(DᵥV,3,3)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DᵥV)[3], codomain(DᵥV)[3])[:,:]
                        end
                    end
                end
            end
        end
    
        for n in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[n,m,ℓ] != 0
                        if n-2 >= 0
                            component(DᵥV,5,1)[:,:] = component(DᵥV,5,1)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ),domain(DᵥV)[1], codomain(DᵥV)[5])[:,:]
                        end
                        if m-2 >= 0
                            component(DᵥV,5,3)[:,:] = component(DᵥV,5,3)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DᵥV)[3], codomain(DᵥV)[5])[:,:]
                        end
                    end
                end
            end
        end
        component(DᵥV,1,1)[:,:] .= 0
        component(DᵥV,1,2)[:,:] .= 0
        component(DᵥV,1,3)[:,:] .= 0
        component(DᵥV,1,4)[:,:] .= 0
        component(DᵥV,1,κ)[:,:] .= 1
     
        
        
    # DᵧV
        DᵧV[:,:] .= 0
    
        codomainᵥ = codomain(DᵧV)[2:5]
        domainᵥ = domain(DᵧV)
        Nᵥ =  order(codomainᵥ)
        Nᵧ =  order(domainᵥ)
    
        Dγ₁V₂temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[2])
        Dγ₃V₂temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[2])
        Dγ₁V₄temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[4])
        Dγ₃V₄temp = zeros(eltype(DᵧV), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[4])
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-3 >= 0 
                            Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  -ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*v₁ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:] 
                        end
                        if n-2 >= 0 && m-2 >= 0
                            Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₃ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:]
                            Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₁ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:]  
                        end
                        if m-3 >=0
                            Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:] 
                        end
                    end
                end
            end
        end
    
        for n in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[n,m,ℓ] != 0
                        if n-3 >= 0 
                            Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*v₁ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:] 
                        end
                        if n-2 >= 0 && m-2 >= 0
                            Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₃ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:]
                            Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*v₁ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:]  
                        end
                        if m-3 >=0
                            Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:] 
                        end
                    end
                end
            end
        end
    
    
        D_even = zeros(eltype(Dγ₁V₂temp) ,2*Nᵧ[1]+1 , Nᵧ[1] +1)
        for j = 0:Nᵧ[1]
            if j == 0
                D_even[1+Nᵧ[1], j+1] = 1
            else
                D_even[1+Nᵧ[1]+j, j+1] = 1
                D_even[1+Nᵧ[1]-j, j+1] = 1
            end
        end
        component(DᵧV, 3 ,1)[:,:] = Dγ₁V₂temp[:,:]*D_even[:,:]
        component(DᵧV, 3 ,3)[:,:] = Dγ₃V₂temp[:,:]*D_even[:,:]
        
        
        D_even = zeros(eltype(Dγ₁V₂temp),2*Nᵧ[3]+1 , Nᵧ[3] +1)
        for j = 0:Nᵧ[3]
            if j == 0
                D_even[1+Nᵧ[3], j+1] = 1
            else
                D_even[1+Nᵧ[3]+j, j+1] = 1
                D_even[1+Nᵧ[3]-j, j+1] = 1
            end
        end
        component(DᵧV, 5 ,1)[:,:] = Dγ₁V₄temp[:,:]*D_even[:,:]
        component(DᵧV, 5 ,3)[:,:] = Dγ₃V₄temp[:,:]*D_even[:,:]   

        #DλC
        component( DλC , 1)[:] .= 0
        component( DλC , 2)[:] = component(c,1)[:] 
        component( DλC , 3)[:] = component(c,2)[:] 
        component( DλC , 4)[:] = component(c,3)[:] 
        component( DλC , 5)[:] = component(c,4)[:] 
        
        
        #DᵥC
        DᵥC[:,:] .= 0
        component( DᵥC , 2 , 1)[:,:] = c₀*Matrix(I, 2*order(codomain(DᵥC)[2])+1, 2*order(domain(DᵥC)[1])+1)
        component( DᵥC , 3 , 2)[:,:] = c₀*Matrix(I, 2*order(codomain(DᵥC)[3])+1, 2*order(domain(DᵥC)[2])+1)
        component( DᵥC , 4 , 3)[:,:] = c₀*Matrix(I, 2*order(codomain(DᵥC)[4])+1, 2*order(domain(DᵥC)[3])+1)
        component( DᵥC , 5 , 4)[:,:] = c₀*Matrix(I, 2*order(codomain(DᵥC)[5])+1, 2*order(domain(DᵥC)[4])+1)

        #DωC
        DωC[:] .= 0
        component(DωC,1)[1] = 0
        component(DωC,2)[:] = project(-c₂,space(DωC)[2])[:]
        component(DωC,3)[:] .= 0
        component(DωC,4)[:] = project(-c₄,space(DωC)[4])[:]
        component(DωC,5)[:] .= 0
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-2 >= 0 && m-2>= 0
                            component(DωC,3)[:] = component(DωC,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃ ,space(DωC)[3])[:]
                        elseif n-2 >= 0
                            component(DωC,3)[:] = component(DωC,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁ ,space(DωC)[3])[:]
                        elseif m-2 >= 0
                            component(DωC,3)[:] = component(DωC,3)[:] - α^(ℓ-1)*M₁[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃ ,space(DωC)[3])[:]
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
                            component(DωC,5)[:] = component(DωC,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁ + (m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃  ,space(DωC)[5])[:]
                        elseif n-2 >= 0
                            component(DωC,5)[:] = component(DωC,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((n-1)*γ₁^(n-2)*γ₃^(m-1)*c₁  ,space(DωC)[5])[:]
                        elseif m-2 >= 0
                            component(DωC,5)[:] = component(DωC,5)[:] - α^(ℓ-1)*M₂[n,m,ℓ] * project((m-1)*γ₁^(n-1)*γ₃^(m-2)*c₃ ,space(DωC)[5])[:]
                        end
                    end
                end
            end
        end
    
        
        #DᵧC
        DᵧC[:,:] .= 0
    
        codomainᵥ = codomain(DᵧC)[2:5]
        domainᵥ = domain(DᵧC)
        Nᵥ =  order(codomainᵥ)
        Nᵧ =  order(domainᵥ)
    
        Dγ₁V₂temp = zeros(eltype(DᵧC), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[2])
        Dγ₃V₂temp = zeros(eltype(DᵧC), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[2])
        Dγ₁V₄temp = zeros(eltype(DᵧC), Fourier(order(domainᵥ[1]), interval(1.0) ),codomainᵥ[4])
        Dγ₃V₄temp = zeros(eltype(DᵧC), Fourier(order(domainᵥ[3]), interval(1.0) ),codomainᵥ[4])
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-3 >= 0 
                            Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  -ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*c₁ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:] 
                        end
                        if n-2 >= 0 && m-2 >= 0
                            Dγ₁V₂temp[:,:] = Dγ₁V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*c₃ ) ,domain(Dγ₁V₂temp), codomain(Dγ₁V₂temp))[:,:]
                            Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*c₁ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:]  
                        end
                        if m-3 >=0
                            Dγ₃V₂temp[:,:] = Dγ₃V₂temp[:,:]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*c₃ ) ,domain(Dγ₃V₂temp), codomain(Dγ₃V₂temp))[:,:] 
                        end
                    end
                end
            end
        end
    
        for n in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[n,m,ℓ] != 0
                        if n-3 >= 0 
                            Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (n-1)*(n-2)*γ₁^(n-3)*γ₃^(m-1)*c₁ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:] 
                        end
                        if n-2 >= 0 && m-2 >= 0
                            Dγ₁V₄temp[:,:] = Dγ₁V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*c₃ ) ,domain(Dγ₁V₄temp), codomain(Dγ₁V₄temp))[:,:]
                            Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(n-1)*γ₁^(n-2)*γ₃^(m-2)*c₁ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:]  
                        end
                        if m-3 >=0
                            Dγ₃V₄temp[:,:] = Dγ₃V₄temp[:,:]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] * project( Multiplication( (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*c₃ ) ,domain(Dγ₃V₄temp), codomain(Dγ₃V₄temp))[:,:] 
                        end
                    end
                end
            end
        end
    
    
        D_even = zeros(eltype(Dγ₁V₂temp) ,2*Nᵧ[1]+1 , Nᵧ[1] +1)
        for j = 0:Nᵧ[1]
            if j == 0
                D_even[1+Nᵧ[1], j+1] = 1
            else
                D_even[1+Nᵧ[1]+j, j+1] = 1
                D_even[1+Nᵧ[1]-j, j+1] = 1
            end
        end
        component(DᵧC, 3 ,1)[:,:] = Dγ₁V₂temp[:,:]*D_even[:,:]
        component(DᵧC, 3 ,3)[:,:] = Dγ₃V₂temp[:,:]*D_even[:,:]
        
        
        D_even = zeros(eltype(Dγ₁V₂temp),2*Nᵧ[3]+1 , Nᵧ[3] +1)
        for j = 0:Nᵧ[3]
            if j == 0
                D_even[1+Nᵧ[3], j+1] = 1
            else
                D_even[1+Nᵧ[3]+j, j+1] = 1
                D_even[1+Nᵧ[3]-j, j+1] = 1
            end
        end
        component(DᵧC, 5 ,1)[:,:] = Dγ₁V₄temp[:,:]*D_even[:,:]
        component(DᵧC, 5 ,3)[:,:] = Dγ₃V₄temp[:,:]*D_even[:,:]   
        
        #Dc₀C 
        component( Dc₀C , 1)[:] .= 0
        component( Dc₀C , 2)[:] = component(v,1)[:] 
        component( Dc₀C , 3)[:] = component(v,2)[:] 
        component( Dc₀C , 4)[:] = component(v,3)[:] 
        component( Dc₀C , 5)[:] = component(v,4)[:] 
        
        #DcC
        DcC[:,:] .= 0
        component(DcC,2,1)[:,:] = Matrix(project(Derivative(1),domain(DcC)[1],codomain(DcC)[2], eltype(DcC) )[:,:]) + project(Multiplication(λ*(c₁^0)),domain(DcC)[1],codomain(DcC)[2])[:,:] 
        component(DcC,2,2)[:,:] = -ω*Matrix(I, 2*order(codomain(DcC)[2])+1, 2*order(domain(DcC)[2])+1)
        component(DcC,2,3)[:,:] = zeros(2*order(codomain(DcC)[2])+1,2*order(domain(DcC)[3])+1)
        component(DcC,2,4)[:,:] = zeros(2*order(codomain(DcC)[2])+1,2*order(domain(DcC)[4])+1)
    
        component(DcC,3,1)[:,:] = zeros(2*order(codomain(DcC)[3])+1,2*order(domain(DcC)[1])+1)
        component(DcC,3,2)[:,:] = Matrix(project(Derivative(1),domain(DcC)[2],codomain(DcC)[3], eltype(DcC))[:,:]) + project(Multiplication(λ*(c₂^0)),domain(DcC)[2],codomain(DcC)[3])[:,:] 
        component(DcC,3,3)[:,:] = zeros(2*order(codomain(DcC)[3])+1,2*order(domain(DcC)[3])+1)
        component(DcC,3,4)[:,:] = zeros(2*order(codomain(DcC)[3])+1,2*order(domain(DcC)[4])+1)
    
        component(DcC,4,1)[:,:] = zeros(2*order(codomain(DcC)[4])+1,2*order(domain(DcC)[1])+1)
        component(DcC,4,2)[:,:] = zeros(2*order(codomain(DcC)[4])+1,2*order(domain(DcC)[2])+1)
        component(DcC,4,3)[:,:] = Matrix(project(Derivative(1),domain(DcC)[3],codomain(DcC)[4], eltype(DcC) )[:,:]) + project(Multiplication(λ*(c₃^0)), domain(DcC)[3],codomain(DcC)[4])[:,:] 
        component(DcC,4,4)[:,:] = -ω*Matrix(I, 2*order(codomain(DcC)[4])+1, 2*order(domain(DcC)[4])+1)
    
        component(DcC,5,1)[:,:] = zeros(2*order(codomain(DcC)[5])+1,2*order(domain(DcC)[1])+1)
        component(DcC,5,2)[:,:] = zeros(2*order(codomain(DcC)[5])+1,2*order(domain(DcC)[2])+1)
        component(DcC,5,3)[:,:] = zeros(2*order(codomain(DcC)[5])+1,2*order(domain(DcC)[3])+1)
        component(DcC,5,4)[:,:] = Matrix(project(Derivative(1),domain(DcC)[4],codomain(DcC)[5], eltype(DcC))[:,:]) + project(Multiplication(λ*(c₄^0)), domain(DcC)[4],codomain(DcC)[5])[:,:] 
    
        for n in axes(M₁,1)
            for m in axes(M₁,2)
                for ℓ in axes(M₁,3)
                    if M₁[n,m,ℓ] != 0
                        if n-2 >= 0
                            component(DcC,3,1)[:,:] = component(DcC,3,1)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ) ,domain(DcC)[1], codomain(DcC)[3])[:,:]
                        end
                        if m-2 >= 0
                            component(DcC,3,3)[:,:] = component(DcC,3,3)[:,:]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DcC)[3], codomain(DcC)[3])[:,:]
                        end
                    end
                end
            end
        end
    
        for n in axes(M₂,1)
            for m in axes(M₂,2)
                for ℓ in axes(M₂,3)
                    if M₂[n,m,ℓ] != 0
                        if n-2 >= 0
                            component(DcC,5,1)[:,:] = component(DcC,5,1)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (n-1)*γ₁^(n-2)*γ₃^(m-1) ),domain(DcC)[1], codomain(DcC)[5])[:,:]
                        end
                        if m-2 >= 0
                            component(DcC,5,3)[:,:] = component(DcC,5,3)[:,:] - α^(ℓ-1)*ω*M₂[n,m,ℓ] * project( Multiplication( (m-1)*γ₁^(n-1)*γ₃^(m-2) ),domain(DcC)[3], codomain(DcC)[5])[:,:]
                        end
                    end
                end
            end
        end
        component(DcC,1,1)[:,:] .= 0
        component(DcC,1,2)[:,:] .= 0
        component(DcC,1,3)[:,:] .= 0
        component(DcC,1,4)[:,:] .= 0
        component(DcC,1,κ)[:,:] .= 1


        return  DλV, DᵥV,DωV,DᵧV, DλC, DᵥC,DωC,DᵧC,Dc₀C , DcC
    end


