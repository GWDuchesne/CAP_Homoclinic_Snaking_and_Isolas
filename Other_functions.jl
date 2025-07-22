# =================================================
# Compilation of others functions used in the proof
# =================================================

# Clear Functions
function clc()
    if Sys.iswindows()
        return read(run(`powershell cls`), String)
    elseif Sys.isunix()
        return read(run(`clear`), String)
    elseif Sys.islinux()
        return read(run(`printf "\033c"`), String)
    end
end

# Extension of γ on exponential Fourier
function γ_extened!(γ)
    Nᵥ = order(γ)

    γ₁ = component(γ,1)[:]
    γ₂ = component(γ,2)[:]
    γ₂ = [0;γ₂]
    γ₃ = component(γ,3)[:]
    γ₄ = component(γ,4)[:]
    γ₄ = [0;γ₄]
    
    γ₁⃰ = zeros(typeof(1im*γ₁[1]),4*Nᵥ[1]+1)
    γ₂⃰ = zeros(typeof(1im*γ₂[1]),4*Nᵥ[2]+1)
    γ₃⃰ = zeros(typeof(1im*γ₃[1]),4*Nᵥ[3]+1)
    γ₄⃰ = zeros(typeof(1im*γ₄[1]),4*Nᵥ[4]+1)
    space_ext = Fourier(2*Nᵥ[1], 1/2) ×Fourier(2*Nᵥ[2], 1/2) ×Fourier(2*Nᵥ[3], 1/2) × Fourier(2*Nᵥ[4], 1/2) 
    

    for i = -2*Nᵥ[1]:2*Nᵥ[1]
        if mod(i,2) == 0
            γ₁⃰[i+2*Nᵥ[1]+1] = γ₁[abs(Int(i/2)) + 1]
        end
    end
    for i = -2*Nᵥ[2]:2*Nᵥ[2]
        if mod(i,2) == 0
            if i == 0
                γ₂⃰[i+2*Nᵥ[2]+1] = γ₂[1]
            else
                γ₂⃰[i+2*Nᵥ[2]+1] = -1im*sign(i)*γ₂[abs(Int(i/2)) + 1]
            end
        end
    end
    for i = -2*Nᵥ[3]:2*Nᵥ[3]
        if mod(i,2) == 0
            γ₃⃰[i+2*Nᵥ[3]+1] = γ₃[abs(Int(i/2)) + 1]
        end
    end
    for i = -2*Nᵥ[4]:2*Nᵥ[4]
        if mod(i,2) == 0
            if i == 0
                γ₄⃰[i+2*Nᵥ[4]+1] = γ₄[1]
            else
                γ₄⃰[i+2*Nᵥ[4]+1] = -1im*sign(i)*γ₄[abs(Int(i/2)) + 1]
            end
        end
    end
    γ⃰ = Sequence( space_ext , [γ₁⃰;γ₂⃰;γ₃⃰; γ₄⃰ ] )
    return γ⃰
end
function γ_extened_orientable!(γ)
    Nᵥ = order(γ)

    γ₁ = component(γ,1)[:]
    γ₂ = component(γ,2)[:]
    zero_type = zeros(typeof(γ₁[1]) ,1,1)[1]
    γ₂ = [zero_type;γ₂]
    γ₃ = component(γ,3)[:]
    γ₄ = component(γ,4)[:]
    γ₄ = [zero_type;γ₄]
    
    γ₁⃰ = Sequence( Fourier(Nᵥ[1], 1), zeros(typeof(1im*γ₁[1]),2*Nᵥ[1]+1))
    γ₂⃰ = Sequence( Fourier(Nᵥ[2], 1), zeros(typeof(1im*γ₂[1]),2*Nᵥ[2]+1))
    γ₃⃰ = Sequence( Fourier(Nᵥ[3], 1), zeros(typeof(1im*γ₃[1]),2*Nᵥ[3]+1))
    γ₄⃰ = Sequence( Fourier(Nᵥ[4], 1), zeros(typeof(1im*γ₄[1]),2*Nᵥ[4]+1))
    space_ext = Fourier(Nᵥ[1], 1) ×Fourier(Nᵥ[2], 1) ×Fourier(Nᵥ[3], 1) × Fourier(Nᵥ[4], 1) 
    
    if typeof(γ₁[1]) == Interval{Float64}
        im_type = interval(1im)
    else
        im_type = 1im
    end
    for i = 0:Nᵥ[1]
        if i == 0
            γ₁⃰[i] = γ₁⃰[i]+γ₁[i+1]
        else
            γ₁⃰[i] = γ₁⃰[i]+γ₁[i+1]
            γ₁⃰[-i] = γ₁⃰[-i]+γ₁[i+1]
        end
    end
    for i = 1:Nᵥ[2]
        γ₂⃰[i] = γ₂⃰[i] - im_type*γ₂[i+1]
        γ₂⃰[-i] = γ₂⃰[-i] + im_type*γ₂[i+1]

    end
    for i = 0:Nᵥ[3]
        if i ==0
            γ₃⃰[i] = γ₃⃰[i] +γ₃[i+1]
        else
            γ₃⃰[i] = γ₃⃰[i] +γ₃[i+1]
            γ₃⃰[-i] = γ₃⃰[-i] + γ₃[i+1]
        end
    end
    for i = 1:Nᵥ[4]
        γ₄⃰[i] = γ₄⃰[i] - im_type*γ₄[i+1]
        γ₄⃰[-i] = γ₄⃰[-i] + im_type*γ₄[i+1]
 
    end

    γ⃰ = Sequence( space_ext , [γ₁⃰[:];γ₂⃰[:];γ₃⃰[:]; γ₄⃰[:]] )
    return γ⃰
end

# Variable Converter
function x2var!(x)
    γ = real(component(x,1:4)) 
    λ = component(x,5)[1]
    v = component(x,6:9)
    w = component(x,10:13)
    L = real(component(x,14)[1])
    θ₀ = real(component(x,15:16))
    σ = real(component(x,17:18))
    a = real(component(x,19:22))
    ω = real(component(x,23)[1])
    eigenpairs = component(x,24)
    p = component(x,25:28)
    return γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p

end

# Remove complex argument for real variable
function x_remove_complex!(x)
    γₙ,λₙ,vₙ,wₙ,Lₙ,θ₀ₙ,σₙ,aₙ,ωₙ,eigenpairsₙ,pₙ = x2var!(x)
    x =  Sequence( space(x),[ γₙ[:]; λₙ ;vₙ[:]; wₙ[:];Lₙ ; θ₀ₙ[:]; σₙ[:]; aₙ[:] ; ωₙ ; eigenpairsₙ[:]; pₙ[:]])
    return x

end

# Space of f
function space_f!(x)
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)
    space_f =  space(Derivative(1)*γ) × ParameterSpace() × space(v) × space(w) ×  ParameterSpace() × ParameterSpace() × ParameterSpace()× ParameterSpace()×  space(a) × ParameterSpace() × ParameterSpace() × ParameterSpace()^10 × space(p)
    return space_f

end

# Space of F
function space_F!(X)
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(component(X,2:29))
    space_F =  ParameterSpace() × space(Derivative(1)*γ) × ParameterSpace() × space(v) × space(w) ×  ParameterSpace() × ParameterSpace() × ParameterSpace()× ParameterSpace()×  space(a) × ParameterSpace() × ParameterSpace() × ParameterSpace()^10 × space(p)
    return space_F

end

# Resize X for padding
function smaller_coeff_in_period(X,Nᵧ,Nᵥ,N_w,Nₐ,Nₚ)

    x = x_remove_complex!(component(X,2:29))
    γ,λ,v,w,L,θ₀,σ,a,ω,eigenpairs,p = x2var!(x)

    space_γ_int_ext = CosFourier( Nᵧ , 1 ) × SinFourier( Nᵧ , 1 ) × CosFourier( Nᵧ, 1 ) × SinFourier(Nᵧ, 1 ) 
    space_v_int_ext = Fourier( Nᵥ , 1 ) × Fourier(Nᵥ , 1 ) × Fourier( Nᵥ , 1 ) × Fourier( Nᵥ , 1 ) 
    space_w_int_ext = (Fourier( N_w[1], 1 ) ⊗ Taylor( N_w[2]  ) ) ×   (Fourier( N_w[1] , 1 ) ⊗ Taylor( N_w[2]  ) ) ×   (Fourier( N_w[1], 1 ) ⊗ Taylor( N_w[2] ) ) ×   (Fourier( N_w[1] , 1 ) ⊗ Taylor( N_w[2]  ) ) 
    space_a_int_ext = Chebyshev(Nₐ) × Chebyshev(Nₐ) × Chebyshev(Nₐ) × Chebyshev(Nₐ) 
    space_p_int_ext = (Taylor( Nₚ[1]) ⊗ Taylor( Nₚ[2])) × (Taylor(Nₚ[1]) ⊗ Taylor(Nₚ[2])) × (Taylor(Nₚ[1]) ⊗ Taylor( Nₚ[2])) × (Taylor( Nₚ[1]) ⊗ Taylor( Nₚ[2])) 

    space_X = ParameterSpace() × space_γ_int_ext × ParameterSpace() × space_v_int_ext × space_w_int_ext × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × ParameterSpace() × space_a_int_ext × ParameterSpace() ×  ParameterSpace()^10 × space_p_int_ext
    return space_X

end

# T operator from Chebyshev 
function T_Chebyshev!( domain_X , codomain_X )
    T = zeros( domain_X, codomain_X )
    for k = 0:order(domain_X)
        for j = 1:order(codomain_X)
            if k == j-1 
                T[j,k] = -1
            end
            if k == j+1 
                T[j,k] = 1
            end
        end
    end
    return T
end

# For syntax 
function turn_exp(strtest, s::AbstractString)
    codes = Dict(collect("1234567890") .=> collect("¹²³⁴⁵⁶⁷⁸⁹⁰"))
    return strtest * map(c -> codes[c], s)
end

function str_system( M₁ )
    strf₁ = ""
    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    if n == 1
                        str_uⁱ = ""
                    elseif n == 2
                        str_uⁱ = "u "
                    else
                        str_uⁱ = turn_exp("u", "$(Int(n-1))")*" "
                    end
                    if m == 1
                        str_vʲ = ""
                    elseif m == 2
                        str_vʲ = "v "
                    else
                        str_vʲ = turn_exp("v", "$(Int(m-1))")*" "
                    end
                    if ℓ == 1
                        str_αℓ = ""
                    elseif ℓ == 2
                        str_αℓ = "α "
                    else
                        str_αℓ = turn_exp("α", "$(Int(ℓ-1))")*" "
                    end
                    if M₁[n,m,ℓ] < 0 || length(strf₁) == 0
                        if M₁[n,m,ℓ] == -1
                            strf₁ = strf₁*String(" - ")*str_αℓ*str_uⁱ*str_vʲ
                        elseif M₁[n,m,ℓ] == 1
                            strf₁ = strf₁*String(" ")*str_αℓ*str_uⁱ*str_vʲ
                        else
                            strf₁ = strf₁*String(" $(M₁[n,m,ℓ]) ")*str_αℓ*str_uⁱ*str_vʲ
                        end
                    else
                        if M₁[n,m,ℓ] == 1
                            strf₁ = strf₁*String(" + ")*str_αℓ*str_uⁱ*str_vʲ
                        else
                            strf₁ = strf₁*String(" + $(M₁[n,m,ℓ]) ")*str_αℓ*str_uⁱ*str_vʲ
                        end
                    end
                end
            end
        end
    end
    return strf₁
end

# Implicit relation of the phase and θ₁ + i θ₂
function expo_to_rad!(a,b)
    if a >= sqrt(2)/2 
        if b > 0
            θ =  atan(b/a)
        else
            θ =  2π - atan(abs(b/a))
        end
    elseif a <= -sqrt(2)/2
        if b > 0
            θ =  π - atan(abs(b/a))
        else
            θ =  π + atan( abs(b/a))
        end
    else
        if b > 0
            if a >0
                θ =  π/2 - atan(abs(a/b))
            else
                θ =  π/2 + atan(abs(a/b))
            end
        else
            if a >0
                θ =  3π/2 + atan(abs(a/b))
            else
                θ =  3π/2 - atan(abs(a/b))
            end
        end
    end
    return θ
end

# Computation of min λ
function λ_min!(λ_cheb)
    N = 100000
    t = range(-1,1,N)
    λ_min = inf.(abs.( λ_cheb(interval(t[1],t[2]))))
    for i = 2:N-1
        λ_min = min(λ_min, inf.(abs.( λ_cheb(interval(t[i],t[i+1])))))
    end
    return λ_min
end

