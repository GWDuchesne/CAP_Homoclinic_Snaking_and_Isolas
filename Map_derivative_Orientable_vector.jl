# DF as vector
function DF_vector!( X, κ, star, lengthᵥ, M₁, M₂ , η )
# Parameters
    α = real(component(X,1)[1])
    x = component(X,2:29)
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
    γ₁, γ₂, γ₃, γ₄ = eachcomponent(γ)
    v₁, v₂, v₃, v₄ = eachcomponent(v)
    w₁, w₂, w₃, w₄ = eachcomponent(w) 
    a₁, a₂, a₃, a₄ = eachcomponent(a) 
    p₁, p₂, p₃, p₄ = eachcomponent(p) 
    σ₁, σ₂ = eachcomponent(σ)
    θ₁, θ₂ = eachcomponent(θ₀)
    γ_ext =  γ_extened_orientable!(γ)
    γ₁_ext = component(γ_ext,1)
    γ₃_ext = component(γ_ext,3)

    DF = Matrix{Any}( undef , 29 , 29 )
    M_index = Vector{Vector}( undef , 29  )
    D₁₀ = Derivative(1,0)
    D₀₁ = Derivative(0,1)

# DPlane
    M_index[1] = []

# DΓ
    # DΓ₁
    M_index[2] = [ 3,24 ]
    DF[ 2 , 3 ] = ω  # Dᵧ₂
    DF[ 2 , 24 ] = γ₂ # D_ω

    # DΓ₂
    M_index[3] = [ 1,2,4,24 ]
    DF[ 3 , 1 ]  = zeros(eltype(γ₁),space(γ₁)) #D_α
    DF[ 3 , 2 ] = zeros(eltype(γ₁),space(γ₁)) # Dᵧ₁
    DF[ 3 , 4 ] = zeros(eltype(γ₁),space(γ₁)) # Dᵧ₃
    DF[ 3 , 24 ] = zeros(eltype(γ₁),space(γ₁)) # D_ω


    # DΓ₃
    M_index[4] = [ 5,24 ]
    DF[ 4 , 5 ] = ω # Dᵧ₄
    DF[ 4 , 24 ] = γ₄ # D_ω

    # DΓ₄
    M_index[5] = [ 1,2 ,4,24 ]
    DF[ 5 , 1 ]  = zeros(eltype(γ₁), space(γ₁)) #D_α
    DF[ 5 , 2 ] = zeros(eltype(γ₁),space(γ₁)) # Dᵧ₁
    DF[ 5 , 4 ] = zeros(eltype(γ₁),space(γ₁)) # Dᵧ₃
    DF[ 5 , 24 ] = zeros(eltype(γ₁),space(γ₁)) # D_ω

# DV
    # DV₀
    M_index[6] = []

    # DV₁
    M_index[7] = [6, 8,24]
    DF[ 7 , 6 ] = v₁    # D_λ
    DF[ 7 , 8 ] = -ω    # Dᵥ₂
    DF[ 7 , 24 ] = -v₂ # D_ω

    # DV₂
    M_index[8] = [1,2,4,6,7,9,24]
    DF[ 8 , 6 ] = v₂    # D_λ
    DF[ 8 , 1 ]  = zeros(eltype(v₁),space(v₁)) # D_α
    DF[ 8 , 2 ]  = zeros(eltype(v₁),space(v₁)) # Dᵧ₁s
    DF[ 8 , 4 ]  = zeros(eltype(v₁),space(v₁)) # Dᵧ₃
    DF[ 8 , 7 ]  = zeros(eltype(v₁),space(v₁)) # Dᵥ₁
    DF[ 8 , 9 ]  = zeros(eltype(v₁),space(v₁)) # Dᵥ₃
    DF[ 8 , 24 ] = zeros(eltype(v₁),space(v₁)) # D_ω

    # DV₃   
    M_index[9] = [6, 10,24]
    DF[ 9 , 6 ] = v₃    # D_λ
    DF[ 9 , 10 ] = -ω  # Dᵥ₄
    DF[ 9 , 24 ] = -v₄ # D_ω

    # DV₄
    M_index[10] = [1,2,4,6,7,9,24]
    DF[ 10 , 6 ] = v₄    # D_λ
    DF[ 10 , 1 ]  = zeros(eltype(v₁), space(v₁)) # D_α
    DF[ 10 , 2 ]  = zeros(eltype(v₁),space(v₁)) # Dᵧ₁
    DF[ 10 , 4 ]  = zeros(eltype(v₁),space(v₁)) # Dᵧ₃
    DF[ 10 , 7 ]  = zeros(eltype(v₁),space(v₁)) # Dᵥ₁
    DF[ 10 , 9 ]  = zeros(eltype(v₁),space(v₁)) # Dᵥ₃
    DF[ 10 , 24 ] = zeros(eltype(v₁),space(v₁)) # D_ω

# DW
    D₀₁ = Derivative(0,1)
    M_FxT = Sequence(Fourier(0,1) ⊗ Taylor(1), [0, 1])

    # DW₁
    M_index[11] = [6,12,24]
    DF_11_24_temp  = M_FxT*(D₀₁*w₁)
    DF_11_24_temp[ (:,0:1) ] .= 0 
    DF[ 11 , 6 ] = DF_11_24_temp # D_λ   

    w₂_temp = copy(w₂)
    w₂_temp[(:, 0:1)] .= 0 
    DF[ 11 , 24 ] = -w₂_temp # D_ω

    DF[ 11 , 12 ] = -ω # D_w₁

    # DW₂
    M_index[12] = [1,6,11,13,24]
    DF[ 12 , 1 ] = zeros( eltype(w₁) , space(w₁)) # D_α


    DF_12_24_temp  = M_FxT*(D₀₁*w₂)
    DF_12_24_temp[ (:,0:1) ] .= 0 
    DF[ 12 , 6 ] = DF_12_24_temp # D_λ   
    DF[ 12 , 24 ] = zeros(eltype(w₁),space(w₁)) # D_ω

    DF[ 12 , 11 ] = zeros( eltype(w₁) , space(w₁)) # D_w₁
    DF[ 12 , 13 ] = zeros( eltype(w₁) , space(w₁)) # D_w₃

    # DW₃
    M_index[13] = [6,14,24]
    DF_13_24_temp  = M_FxT*(D₀₁*w₃)
    DF_13_24_temp[ (:,0:1) ] .= 0 
    DF[ 13 , 6 ] = DF_13_24_temp # D_λ   

    w₄_temp = copy(w₄)
    w₄_temp[(:, 0:1)] .= 0 
    DF[ 13 , 24 ] = -w₄_temp # D_ω

    DF[ 13 , 14 ] = -ω # D_w₄

    # DW₄
    M_index[14] = [1,6,11,13,24]
    DF[ 14 , 1 ] = zeros( eltype(w₁) , space(w₁)) # D_α

    DF_14_24_temp  = M_FxT*(D₀₁*w₄)
    DF_14_24_temp[ (:,0:1) ] .= 0 
    DF[ 14 , 6 ] = DF_14_24_temp # D_λ   
    DF[ 14 , 24 ] = zeros(eltype(w₁),space(w₁)) # D_ω

    DF[ 14 , 11 ] = zeros( eltype(w₁) , space(w₁)) # D_w₁
    DF[ 14 , 13 ] = zeros( eltype(w₁) , space(w₁)) # D_w₃
# DG
    # DH₁
    M_index[15] = []

    # DH₂
    M_index[16] = []

    # DH₃
    M_index[17] = []

    # DH₄
    M_index[18] = []

    for i = 1:4
        DF[14+i,18] = (D₁₀*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) + (D₀₁*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )
        DF[14+i,19] = 1im*(D₁₀*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] ) - 1im*(D₀₁*component(p,i))( σ₁[1] + 1im*σ₂[1], σ₁[1] - 1im*σ₂[1] )    
    end 

    DF_19to22_16 = Vector{Any}(undef, 4 )
    DF_19to22_17 = Vector{Any}(undef, 4 )
    for i = 1:4
        N_w = order(w)[i] 
        wᵢ = component(w,i)
        n = Sequence(space(wᵢ),reshape(repeat(-N_w[1]:N_w[1],1,N_w[2]+1), (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        m = Sequence(space(wᵢ),reshape(repeat(transpose(0:N_w[2]),2*N_w[1]+1,1) , (2*N_w[1]+1)*(N_w[2]+1)  ,1)[:])
        DF_19to22_16[i] = sum(n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)] + abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
        DF_19to22_17[i] = sum(1im*n[(1:N_w[1],:)].*wᵢ[(1:N_w[1],:)].*((θ₁[1] + 1im*θ₂[1]).^(n[(1:N_w[1],:)].-1)).*(-0.95).^m[(1:N_w[1],:)]  - 1im*abs.(n[(-N_w[1]:-1,:)]).*wᵢ[(-N_w[1]:-1,:)].*((θ₁[1] - 1im*θ₂[1]).^(abs.(n[(-N_w[1]:-1,:)]) .-1)).*(-0.95).^m[(-N_w[1]:-1,:)] ) 
   end


    # DG₁
    M_index[19] = [15,16,17,21]
    Ψ_19_15 = 1/2*a₂
    DF_19_16 = zeros( eltype(w₁),  space(a₁)  )
    DF_19_16[0] = DF_19to22_16[1]
    DF[19,16] = DF_19_16
    DF_19_17 = zeros(   eltype(w₁),  space(a₁) )
    DF_19_17[0] = DF_19to22_17[1]
    DF[19,17] = DF_19_17

    Ψ_19_21 = Sequence( Chebyshev(0) , [L/2] )

    # DG₂
    M_index[20] = [1, 15,16,17,20,22]
    DαΨ₂ = zeros( eltype(w₁),  space(a₁) )

    Ψ_20_15 = zeros( eltype(w₁),  space(a₁) )
    DF_20_16 = zeros(  eltype(w₁),  space(a₂)  )
    DF_20_16[0] = DF_19to22_16[2]
    DF[20,16] = DF_20_16
    DF_20_17 = zeros(   eltype(w₁),  space(a₂ ))
    DF_20_17[0] = DF_19to22_17[2]
    DF[20,17] = DF_20_17

    Ψ_20_20 = zeros(   eltype(w₁),  space(a₁) )
    Ψ_20_22 = zeros(   eltype(w₁),  space(a₁) )

    # DG₃
    M_index[21] = [ 15,16,17,23]
    Ψ_21_15 = 1/2*a₄
    DF_21_16 = zeros(   eltype(w₁),  space(a₃ ) )
    DF_21_16[0] = DF_19to22_16[3]
    DF[21,16] = DF_21_16
    DF_21_17 = zeros(   eltype(w₁),  space(a₃) )
    DF_21_17[0] = DF_19to22_17[3]
    DF[21,17] = DF_21_17

    Ψ_21_23  =  Sequence( Chebyshev(0) , [L/2] )

    # DG₄
    M_index[22] = [1, 15,16,17,20,22]
    DαΨ₄ = zeros( eltype(w₁),  space(a₁) )
    Ψ_22_15 = zeros( eltype(w₁),  space(a₁) )
    DF_22_16 = zeros(   eltype(w₁),  space(a₄ ) )
    DF_22_16[0] = DF_19to22_16[4]
    DF[22,16] = DF_22_16
    DF_22_17 = zeros(   eltype(w₁),  space(a₄)  )
    DF_22_17[0] = DF_19to22_17[4]
    DF[22,17] = DF_22_17

    Ψ_22_20 = zeros(   eltype(w₁),  space(a₁) )
    Ψ_22_22 = zeros(   eltype(w₁),  space(a₁) )

    # DQ₁
    M_index[23] = [18,19]
    DF[ 23 , 18 ] = -2*σ₁
    DF[ 23 , 19 ] = -2*σ₂

    # DQ₂
    M_index[24] = [ 16, 17]
    DF[ 24 , 16 ] = -2*θ₁
    DF[ 24 , 17 ] = -2*θ₂

# DE
    M_index[25] = [1]
    DF[ 25 , 25 ] = DE!(zeros( eltype(eigenpairs) , ParameterSpace()^10, ParameterSpace()^10 ), eigenpairs, M₁, M₂,α, star)
    M = zeros(eltype(w₁),4,4)
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
    DF[ 25 , 1 ]  =  Sequence(ParameterSpace()^10, [ 0;- M*component(eigenpairs,2:5)[:]; 0; - M*component(eigenpairs,7:10)[:] ])

# DP
    # DP₁
    M_index[26] = [ 27 ]
    DF[ 26 , 27 ] = -1 # Dₚ₂

    # DP₁
    M_index[27] = [ 1, 26, 28 ]
    DF[27, 1 ]  = zeros(eltype(p₁),space(p₁)) # D_α

    DF[27, 26 ]  = zeros(eltype(p₁),space(p₁)) # Dₚ₁
    DF[27, 28 ]  = zeros(eltype(p₁),space(p₁)) # Dₚ₃
    # DP₁
    M_index[28] = [ 29 ]
    DF[ 28 , 29 ] = -1 # Dₚ₄

    # DP₁
    M_index[29] = [ 1, 26, 28 ]
    DF[29, 1 ]  = zeros(eltype(p₁),space(p₁)) # D_α 

    DF[29, 26 ]  = zeros(eltype(p₁),space(p₁)) # Dₚ₁
    DF[29, 28 ]  = zeros(eltype(p₁),space(p₁)) # Dₚ₃

# All the loops goes here 
    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    DF[ 3 , 24 ] = DF[ 3 , 24 ] - α^(ℓ-1)*M₁[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1))
                    DF[ 12 , 24 ] = DF[ 12 , 24 ] - α^(ℓ-1)*M₁[n,m,ℓ]*w₁^(n-1)*w₃^(m-1)
                    Ψ_20_15 = Ψ_20_15 + 1/2*M₁[n,m,ℓ] * α^(ℓ-1)*a₁^(n-1)*a₃^(m-1)

                    if n-2 >= 0 && m-2>= 0
                        DF[ 8 , 24 ] = DF[ 8 , 24 ] - α^(ℓ-1)*M₁[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                    elseif n >= 2
                        DF[ 8 , 24 ] = DF[ 8 , 24 ] - α^(ℓ-1)*M₁[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ 
                    elseif m >= 2
                        DF[ 8 , 24 ] = DF[ 8 , 24 ] - α^(ℓ-1)*M₁[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                    end
                    if n >= 2
                        DF[ 3 , 2 ] = DF[ 3 , 2 ] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(n-1)*(γ₁^(n-2)*γ₃^(m-1))
                        DF[ 8 , 7 ] = DF[ 8 , 7 ]  - α^(ℓ-1)*ω*M₁[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) 
                        DF[ 12 , 11 ]  = DF[ 12 , 11 ]  - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(n-1)*w₁^(n-2)*w₃^(m-1)
                        DF[27, 26 ]  = DF[27, 26 ]  - α^(ℓ-1)*M₁[n,m,ℓ]*(n-1)*p₁^(n-2)*p₃^(m-1)
                        Ψ_20_20 = Ψ_20_20 + α^(ℓ-1)*L/2*M₁[n,m,ℓ] * (n-1)*a₁^(n-2)*a₃^(m-1)

                    end
                    if m >= 2
                        DF[ 3 , 4 ] = DF[ 3 , 4 ] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(m-1)*(γ₁^(n-1)*γ₃^(m-2))
                        DF[ 8 , 9 ] = DF[ 8 , 9 ]   - α^(ℓ-1)*ω*M₁[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) 
                        DF[ 12 , 13 ] = DF[ 12 , 13 ] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(m-1)*w₁^(n-1)*w₃^(m-2)
                        DF[27, 28 ]  = DF[27, 28 ]  - α^(ℓ-1)*M₁[n,m,ℓ]*(m-1)*p₁^(n-1)*p₃^(m-2)
                        Ψ_20_22 = Ψ_20_22 + α^(ℓ-1)*L/2*M₁[n,m,ℓ] * (m-1)*a₁^(n-1)*a₃^(m-2)

                    end
                    if n >= 3 
                        DF[ 8 , 2 ] = DF[ 8 , 2 ]  -ω*α^(ℓ-1)*M₁[n,m,ℓ] * (n-1)*(n-2)*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ 
                    end
                    if n >= 2 && m >= 2
                        DF[ 8 , 2 ] = DF[ 8 , 2 ]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] *  (m-1)*(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ 
                        DF[ 8 , 4 ] = DF[ 8 , 4 ]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] *  (m-1)*(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ 
                    end
                    if m >= 3
                        DF[ 8 , 4 ] = DF[ 8 , 4 ]  - ω*α^(ℓ-1)*M₁[n,m,ℓ] *  (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ 
                    end
                    if ℓ >= 2
                        DF[ 3 , 1 ]  = DF[ 3 , 1 ]  - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * γ₁^(n-1)*γ₃^(m-1)
                        DF[ 12 , 1 ] = DF[ 12 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * w₁^(n-1)*w₃^(m-1)
                        DαΨ₂ = DαΨ₂ + L/2*(ℓ-1)*α^(ℓ-2)*M₁[n,m,ℓ] * a₁^(n-1)*a₃^(m-1)
                        DF[27, 1 ]  = DF[27, 1 ]  - (ℓ-1)*α^(ℓ-2)*M₁[n,m,ℓ] * p₁^(n-1)*p₃^(m-1)

                        if n-2 >= 0 && m-2>= 0
                            DF[ 8 , 1 ] = DF[ 8 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                        elseif n-2 >= 0
                            DF[ 8 , 1 ] = DF[ 8 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ 
                        elseif m-2 >= 0
                            DF[ 8 , 1 ] = DF[ 8 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₁[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                        end
                    end

                end
            end
        end
    end

    for n  in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    DF[ 5 , 24 ] = DF[ 5 , 24 ] - α^(ℓ-1)*M₂[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1))
                    DF[ 14 , 24 ] = DF[ 14 , 24 ] - α^(ℓ-1)*M₂[n,m,ℓ]*w₁^(n-1)*w₃^(m-1)
                    Ψ_22_15 = Ψ_22_15 + 1/2*M₂[n,m,ℓ] * α^(ℓ-1)*a₁^(n-1)*a₃^(m-1)

                    if n-2 >= 0 && m-2>= 0
                        DF[ 10 , 24 ] = DF[ 10 , 24 ] - α^(ℓ-1)*M₂[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                    elseif n >= 2
                        DF[ 10 , 24 ] = DF[ 10 , 24 ] - α^(ℓ-1)*M₂[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ 
                    elseif m >= 2
                        DF[ 10 , 24 ] = DF[ 10 , 24 ] - α^(ℓ-1)*M₂[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                    end
                    if n >= 2
                        DF[ 5 , 2 ] = DF[ 5 , 2 ] - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(n-1)*(γ₁^(n-2)*γ₃^(m-1))       
                        DF[ 10 , 7 ] = DF[ 10 , 7 ]  - α^(ℓ-1)*ω*M₂[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) 
                        DF[ 14 , 11 ]  = DF[ 14 , 11 ]  - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(n-1)*w₁^(n-2)*w₃^(m-1)
                        DF[29, 26 ]  = DF[29, 26 ]  - α^(ℓ-1)*M₂[n,m,ℓ]*(n-1)*p₁^(n-2)*p₃^(m-1)
                        Ψ_22_20 = Ψ_22_20 + α^(ℓ-1)*L/2*M₂[n,m,ℓ] * (n-1)*a₁^(n-2)*a₃^(m-1)

                    end
                    if m >= 2
                        DF[ 5 , 4 ] = DF[ 5 , 4 ] - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(m-1)*(γ₁^(n-1)*γ₃^(m-2))
                        DF[ 10 , 9 ] = DF[ 10 , 9 ]   - α^(ℓ-1)*ω*M₂[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) 
                        DF[ 14 , 13 ] = DF[ 14 , 13 ] - α^(ℓ-1)*ω*M₂[n,m,ℓ]*(m-1)*w₁^(n-1)*w₃^(m-2)
                        DF[29, 28 ]  = DF[29, 28 ]  - α^(ℓ-1)*M₂[n,m,ℓ]*(m-1)*p₁^(n-1)*p₃^(m-2)
                        Ψ_22_22 = Ψ_22_22 + α^(ℓ-1)*L/2*M₂[n,m,ℓ] * (m-1)*a₁^(n-1)*a₃^(m-2)

                    end
                    if n >= 3 
                        DF[ 10 , 2 ] = DF[ 10 , 2 ]  -ω*α^(ℓ-1)*M₂[n,m,ℓ] * (n-1)*(n-2)*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ 
                    end
                    if n >= 2 && m >= 2
                        DF[ 10 , 2 ] = DF[ 10 , 2 ]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] *  (m-1)*(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ 
                        DF[ 10 , 4 ] = DF[ 10 , 4 ]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] *  (m-1)*(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ 
                    end
                    if m >= 3
                        DF[ 10 , 4 ] = DF[ 10 , 4 ]  - ω*α^(ℓ-1)*M₂[n,m,ℓ] *  (m-1)*(m-2)*γ₁^(n-1)*γ₃^(m-3)*v₃ 
                    end
                    if ℓ >= 2
                        DF[ 5 , 1 ]  = DF[ 5 , 1 ]  - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * γ₁^(n-1)*γ₃^(m-1)
                        DF[ 14 , 1 ] = DF[ 14 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * w₁^(n-1)*w₃^(m-1)
                        DαΨ₄ = DαΨ₄ + L/2*(ℓ-1)*α^(ℓ-2)*M₂[n,m,ℓ] * a₁^(n-1)*a₃^(m-1)
                        DF[29, 1 ]  = DF[29, 1 ]  - (ℓ-1)*α^(ℓ-2)*M₂[n,m,ℓ] * p₁^(n-1)*p₃^(m-1)

                        if n-2 >= 0 && m-2>= 0
                            DF[ 10 , 1 ] = DF[ 10 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                        elseif n-2 >= 0
                            DF[ 10 , 1 ] = DF[ 10 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * (n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ 
                        elseif m-2 >= 0
                            DF[ 10 , 1 ] = DF[ 10 , 1 ] - (ℓ-1)*α^(ℓ-2)*ω*M₂[n,m,ℓ] * (m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ 
                        end
                    end
                end
            end
        end
    end

    DF_14_1 = copy(DF[ 14 , 1 ])
    DF_14_1[(:,0:1)] .= 0 
    DF[ 14 , 1 ] = DF_14_1

    DF_12_1 = copy(DF[ 12 , 1 ])
    DF_12_1[(:,0:1)] .= 0 
    DF[ 12 , 1 ] = DF_12_1

    DF_14_24 = copy(DF[ 14 , 24 ])
    DF_14_24[(:,0:1)] .= 0 
    DF[ 14 , 24 ] = DF_14_24

    DF_12_24 = copy(DF[ 12 , 24 ])
    DF_12_24[(:,0:1)] .= 0 
    DF[ 12 , 24 ] = DF_12_24

    DF_27_1 = copy(DF[ 27 , 1 ])
    DF_27_1[(0,0)] .= 0 
    DF_27_1[(1,0)] .= 0 
    DF_27_1[(0,2)] .= 0 
    DF[ 27 , 1 ] = DF_27_1

    DF_29_1 = copy(DF[ 29 , 1 ])
    DF_29_1[(0,0)] .= 0 
    DF_29_1[(1,0)] .= 0 
    DF_29_1[(0,1)] .= 0 
    DF[ 29 , 1 ] = DF_29_1

    Ψ_19_15 = project( Ψ_19_15 , Chebyshev(order(Ψ_19_15)+2) ) 
    Ψ_20_15 = project( Ψ_20_15 , Chebyshev(order(Ψ_20_15)+2) ) 
    Ψ_21_15 = project( Ψ_21_15 , Chebyshev(order(Ψ_21_15)+2) ) 
    Ψ_22_15 = project( Ψ_22_15 , Chebyshev(order(Ψ_22_15)+2) ) 

    DF_19_15 = zeros( eltype(Ψ_19_15) ,  Chebyshev(order(Ψ_19_15) -1) )
    DF_20_15 = zeros( eltype(Ψ_20_15) ,  Chebyshev(order(Ψ_20_15) -1) )
    DF_21_15 = zeros( eltype(Ψ_21_15) ,  Chebyshev(order(Ψ_21_15) -1) )
    DF_22_15 = zeros( eltype(Ψ_22_15) ,  Chebyshev(order(Ψ_22_15) -1) )

    DF_19_15[1:end] = DF_19_15[1:end] + Ψ_19_15[2:end]  - Ψ_19_15[0:end-2] 
    DF_20_15[1:end] = DF_20_15[1:end] + Ψ_20_15[2:end]  - Ψ_20_15[0:end-2] 
    DF_21_15[1:end] = DF_21_15[1:end] + Ψ_21_15[2:end]  - Ψ_21_15[0:end-2] 
    DF_22_15[1:end] = DF_22_15[1:end] + Ψ_22_15[2:end]  - Ψ_22_15[0:end-2] 

    DF[19,15] = DF_19_15
    DF[20,15] = DF_20_15
    DF[21,15] = DF_21_15
    DF[22,15] = DF_22_15

    DαΨ₂ = project( DαΨ₂ , Chebyshev(order(DαΨ₂)+2) ) 
    DαΨ₄ = project( DαΨ₄ , Chebyshev(order(DαΨ₄)+2) ) 

    DF_20_1 = zeros( eltype(DαΨ₂) ,  Chebyshev(order(DαΨ₂) -1) )
    DF_22_1 = zeros( eltype(DαΨ₄) ,  Chebyshev(order(DαΨ₄) -1) )

    DF_20_1[1:end] = DF_20_1[1:end] + DαΨ₂[2:end]  - DαΨ₂[0:end-2] 
    DF_22_1[1:end] = DF_22_1[1:end] + DαΨ₄[2:end]  - DαΨ₄[0:end-2] 

    DF[20,1] = DF_20_1
    DF[22,1] = DF_22_1
    
    #=
    Ψ_19_21 = project( Ψ_19_21 , Chebyshev(order(Ψ_19_21)+2) ) 
    Ψ_20_20 = project( Ψ_20_20 , Chebyshev(order(Ψ_20_20)+2) ) 
    Ψ_20_22 = project( Ψ_20_22 , Chebyshev(order(Ψ_20_22)+2) ) 
    Ψ_21_23 = project( Ψ_21_23 , Chebyshev(order(Ψ_21_23)+2) ) 
    Ψ_22_20 = project( Ψ_22_20 , Chebyshev(order(Ψ_22_20)+2) ) 
    Ψ_22_22 = project( Ψ_22_22 , Chebyshev(order(Ψ_22_22)+2) ) 

    DF_19_21 = zeros( eltype(Ψ_19_21) ,  Chebyshev(order(Ψ_19_21) -1) )
    DF_20_20 = zeros( eltype(Ψ_20_20) ,  Chebyshev(order(Ψ_20_20) -1) )
    DF_20_22 = zeros( eltype(Ψ_20_22) ,  Chebyshev(order(Ψ_20_22) -1) )
    DF_21_23 = zeros( eltype(Ψ_21_23) ,  Chebyshev(order(Ψ_21_23) -1) )
    DF_22_20 = zeros( eltype(Ψ_22_20) ,  Chebyshev(order(Ψ_22_20) -1) )
    DF_22_22 = zeros( eltype(Ψ_22_22) ,  Chebyshev(order(Ψ_22_22) -1) )

    DF_19_21[1:end] = DF_19_21[1:end] + Ψ_19_21[2:end]  - Ψ_19_21[0:end-2] 
    DF_20_20[1:end] = DF_20_20[1:end] + Ψ_20_20[2:end]  - Ψ_20_20[0:end-2] 
    DF_20_22[1:end] = DF_20_22[1:end] + Ψ_20_22[2:end]  - Ψ_20_22[0:end-2] 
    DF_21_23[1:end] = DF_21_23[1:end] + Ψ_21_23[2:end]  - Ψ_21_23[0:end-2] 
    DF_22_20[1:end] = DF_22_20[1:end] + Ψ_22_20[2:end]  - Ψ_22_20[0:end-2] 
    DF_22_22[1:end] = DF_22_22[1:end] + Ψ_22_22[2:end]  - Ψ_22_22[0:end-2] 
    =#
    DF[19,21] = Ψ_19_21
    DF[20,20] = Ψ_20_20
    DF[20,22] = Ψ_20_22
    DF[21,23] = Ψ_21_23
    DF[22,20] = Ψ_22_20
    DF[22,22] = Ψ_22_22



return DF,M_index

end
#=
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

    T = zeros( codomain(Da₁Ψ₂), codomain(Da₁Ψ₂) )
    for k = 0:order(codomain(Da₁Ψ₂))
        for j = 1:order(codomain(Da₁Ψ₂))
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
    #=

    space_domain₁ = Chebyshev(order(domain(DₐG)[1])+1) 
    space_domain₃ = Chebyshev(order(domain(DₐG)[3])+1) 
    for i = 5:8
    component(DₐG,i,i-4)[ : ,:] = 2*diagm(0:max(order(domain(DₐG)[i-4]),order(codomain(DₐG)[i]) ) )[ 1:order(codomain(component(DₐG,i,i-4)))+1 ,1:order(domain(component(DₐG,i,i-4)))+1]
    end
    a₁_EXT = project( a₁ , Chebyshev(Nₐ[1]+1) )
    a₃_EXT = project( a₃ , Chebyshev(Nₐ[3]+1) )

    #
    # First the derivative by a
    # 

    Da₂Ψ₁ = L/2*project( Multiplication(a₁_EXT^0 ),space(a₁_EXT),space(a₁_EXT) )
    Da₁Ψ₂ = zeros( eltype(w), space_domain₁,space(a₁_EXT) )[:,:]
    Da₃Ψ₂ = zeros( eltype(w), space_domain₃,space(a₁_EXT) )[:,:]
    Da₄Ψ₃ = L/2*project( Multiplication(a₁_EXT^0 ),space(a₁_EXT),space(a₁_EXT) )
    Da₁Ψ₄ = zeros( eltype(w), space_domain₁,space(a₁_EXT) )[:,:]
    Da₃Ψ₄ = zeros( eltype(w), space_domain₃,space(a₁_EXT) )[:,:]

    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0
                    if i >= 2
                    Da₁Ψ₂ = Da₁Ψ₂ + α^(ℓ-1)*L/2*M₁[i,j,ℓ] * project(  Multiplication((i-1)*a₁_EXT^(i-2)*a₃_EXT^(j-1))   ,space_domain₁,space(a₁_EXT))[:,:]
                    end
                    if j>= 2
                    Da₃Ψ₂ = Da₃Ψ₂ + α^(ℓ-1)*L/2*M₁[i,j,ℓ] * project(  Multiplication((j-1)*a₁_EXT^(i-1)*a₃_EXT^(j-2))   ,space_domain₃,space(a₁_EXT))[:,:]
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
                    Da₁Ψ₄ = Da₁Ψ₄ + α^(ℓ-1)*L/2*M₂[i,j,ℓ] * project(  Multiplication((i-1)*a₁_EXT^(i-2)*a₃_EXT^(j-1))   ,space_domain₁,space(a₁_EXT))[:,:]
                    end
                    if j>= 2
                    Da₃Ψ₄ = Da₃Ψ₄ + α^(ℓ-1)*L/2*M₂[i,j,ℓ] * project(  Multiplication((j-1)*a₁_EXT^(i-1)*a₃_EXT^(j-2))   ,space_domain₃,space(a₁_EXT))[:,:]
                    end
                end
            end
        end
    end
    =#
     #=
    T = zeros( domain(Da₂Ψ₁), codomain(Da₂Ψ₁) )
    for k = 0:order(domain(Da₂Ψ₁))
        for j = 1:order(codomain(Da₂Ψ₁))
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end

   
    T = Matrix( Tridiagonal( -ones(length(a₁),1)[:] ,zeros(length(a₁_EXT),1)[:] , ones(length(a₁),1)[:] )   )
    T[1,:] .= 0

    component(DₐG, 5,2 )[:,:] = (T*Da₂Ψ₁)[1:end-1,1:end-1] 
    component(DₐG, 6,1 )[:,:] = (T*Da₁Ψ₂)[1:end-1,1:end-1] 
    component(DₐG, 6,3 )[:,:] = (T*Da₃Ψ₂)[1:end-1,1:end-1] 
    component(DₐG, 7,4 )[:,:] = (T*Da₄Ψ₃)[1:end-1,1:end-1]
    component(DₐG, 8,1 )[:,:] = (T*Da₁Ψ₄)[1:end-1,1:end-1] 
    component(DₐG, 8,3 )[:,:] = (T*Da₃Ψ₄)[1:end-1,1:end-1] 
    
=#
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
    M = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    for  ℓ in axes(M₁,3)
        M = M + [0 0 0 0; α^(ℓ-1)*M₁_temp[2,1,ℓ] 0 α^(ℓ-1)*M₁_temp[1,2,ℓ] 0; 0 0 0 0; 0 0 0 0]
    end
    for  ℓ in axes(M₂,3)
        M = M + [0 0 0 0; 0 0 0 0; 0 0 0 0; α^(ℓ-1)*M₂_temp[2,1,ℓ] 0 α^(ℓ-1)*M₂_temp[1,2,ℓ] 0]
    end

    δₖ₁ = 0
    δₖ₂ = 0
    δₖ₃ = 0
    δₖ₄ = 0

    if star == 1
        δₖ₁ = 1
    elseif  star == 2
        δₖ₂ = 1
    elseif  star == 3
        δₖ₃ = 1
    elseif  star == 4
        δₖ₄ = 1
    end


    DE[1,:] = [ 0 δₖ₁ δₖ₂ δₖ₃ δₖ₄ 0 0 0 0 0]
    DE[2:5,:] = [ξ₁ 1.0*Matrix(I,4,4)*λ₁-M zeros(4,5) ]   
    DE[6,:] = [ 0 0 0 0 0 0 δₖ₁ δₖ₂ δₖ₃ δₖ₄]
    DE[7:10,:] = [zeros(4,5) ξ₂ 1.0*Matrix(I,4,4)*λ₂-M ]  
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
                        component(DᵧΓ,2,1)[:,:] = component(DᵧΓ,2,1)[:,:] - α^(ℓ-1)*ω*M₁[n,m,ℓ]*(n-1)*project( Multiplication(γ₁^(n-2)*γ₃^(m-1)),domain(DᵧΓ)[3],codomain(DᵧΓ)[1])[:,:]
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


    for i in axes(M₁,1)
        for j in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[i,j,ℓ] != 0 && ℓ > 1
                    component(DαF,2)[:,:]  = component(DαF,2)[:,:]  - (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project(γ₁^(i-1)*γ₃^(j-1),codomain(component(DαF,2 )))[:]
                    component(DαF,11)[:,:] = component(DαF,11)[:,:]- (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project(w₁^(i-1)*w₃^(j-1),codomain(component(DαF,11)))[:]
                    component(DαF,26)[:,:] = component(DαF,26)[:,:] - (ℓ-1)*α^(ℓ-2)*M₁[i,j,ℓ] * project(p₁^(i-1)*p₃^(j-1),codomain(component(DαF,26)))[:]
                    DαΨ₂ = DαΨ₂ + L/2*(ℓ-1)*α^(ℓ-2)*M₁[i,j,ℓ] * project(a₁^(i-1)*a₃^(j-1),space(DαΨ₂))
                    if i-2 >= 0 && j-2>= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project((i-1)*γ₁_ext^(i-2)*γ₃_ext^(j-1)*v₁ + (j-1)*γ₁_ext^(i-1)*γ₃_ext^(j-2)*v₃ ,codomain(component(DαF,7)))[:]
                    elseif i-2 >= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project((i-1)*γ₁_ext^(i-2)*γ₃_ext^(j-1)*v₁ ,codomain(component(DαF,7)))[:]
                    elseif j-2 >= 0
                        component(DαF,7)[:,:] = component(DαF,7)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₁[i,j,ℓ] * project((j-1)*γ₁_ext^(i-1)*γ₃_ext^(j-2)*v₃ ,codomain(component(DαF,7)))[:]
                    end
                end
            end
        end
    end

    for i  in axes(M₂,1)
        for j in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[i,j,ℓ] != 0 && ℓ > 1
                    component(DαF,4 )[:,:] = component(DαF,4 )[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project(γ₁^(i-1)*γ₃^(j-1),codomain(component(DαF,4 )))[:]
                    component(DαF,13)[:,:] = component(DαF,13)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project(w₁^(i-1)*w₃^(j-1),codomain(component(DαF,13)))[:]
                    component(DαF,28)[:,:] = component(DαF,28)[:,:] - (ℓ-1)*α^(ℓ-2)*M₂[i,j,ℓ] * project(p₁^(i-1)*p₃^(j-1),codomain(component(DαF,28)))[:]
                    DαΨ₄ = DαΨ₄ + L/2*(ℓ-1)*α^(ℓ-2)*M₂[i,j,ℓ] * project(a₁^(i-1)*a₃^(j-1),space(DαΨ₄))
                    if i-2 >= 0 && j-2>= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project((i-1)*γ₁_ext^(i-2)*γ₃_ext^(j-1)*v₁ + (j-1)*γ₁_ext^(i-1)*γ₃_ext^(j-2)*v₃ ,codomain(component(DαF,9)))[:]
                    elseif i-2 >= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project((i-1)*γ₁_ext^(i-2)*γ₃_ext^(j-1)*v₁ ,codomain(component(DαF,9)))[:]
                    elseif j-2 >= 0
                        component(DαF,9)[:,:] = component(DαF,9)[:,:] - (ℓ-1)*α^(ℓ-2)*ω*M₂[i,j,ℓ] * project((j-1)*γ₁_ext^(i-1)*γ₃_ext^(j-2)*v₃ ,codomain(component(DαF,9)))[:]
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
=#