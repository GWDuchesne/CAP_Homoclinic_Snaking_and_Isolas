# ===============================================================
# Collection of functions used in the computations of the bounds
# ===============================================================

function A_fft!( A_fft , space_X, space_F , X_fft, κ, star, lengthᵥ, M₁, M₂ , η_fft ,N_fft)
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N_fft]
    for j = 1:N_fft
        DF_j = inv(DF_all!(DF[j], X_fft[j], κ, star, lengthᵥ, M₁, M₂ , η_fft[j] ))
        for i = 1:29
            if j == 1
                A_fft[i] = Vector{Matrix}(undef, N_fft)
            end
            A_fft[i][j] = component( DF_j , : , i)[:,:]
        end
    end
    return A_fft
end


function A_grid!( A_fft , space_X, space_F , X_fft, κ, star, lengthᵥ, M₁, M₂ , η_fft ,N)
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N+1]
    for j = 1:N+1
        DF_j = inv(DF_all!(DF[j], X_fft[j], κ, star, lengthᵥ, M₁, M₂ , η_fft[j] ))
        for i = 1:29
            if j == 1
                A_fft[i] = Vector{Matrix}(undef, N+1)
            end
            A_fft[i][j] = component( DF_j , : , i)[:,:]
        end
    end
    return A_fft
end




function A_cheb!( A_cheb::Vector,  A_fft::Vector{<:Vector} , N::Int64  )
    iter = ProgressBar(1:29)
    set_description(iter, "    - ")
    for i in iter
       A_cheb[i] =  map(A -> interval.(A), grid2cheb(A_fft[i],N))
    end
    return A_cheb
end
function A_cheb_mid!( A_cheb,  A_fft , N  )
    for i = 1:29
        A_cheb[i] =   grid2cheb(A_fft[i],N)
    end
    return A_cheb
end



function A_Taylor!( A_Taylor,  A_grid , N  )
    for i = 1:29
        A_Taylor[i] =  map(A -> interval.(A), grid2taylor(A_grid[i]))
    end
    return A_Taylor
end


function dN!( N , M₁ , M₂ , order_p , order_w)
    M₁_index =  M₁.*0
    M₁_index_α₁₀ =  0
    M₁_index_α₀₁ =  0
    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    M₁_index[n,m,ℓ] = (n-1) + (m-1) + (ℓ-1)
                    if n == 2 && m ==1 && ℓ >= 2
                        M₁_index_α₁₀ = M₁_index_α₁₀ +1 
                    end
                    if m == 2 && n ==1 && ℓ >= 2
                        M₁_index_α₀₁ = M₁_index_α₀₁ +1 
                    end
                end
            end
        end
    end    
    M₂_index =  M₂.*0
    M₂_index_α₁₀ =  0
    M₂_index_α₀₁ =  0
    for n in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    M₂_index[n,m,ℓ] = (n-1) + (m-1) + (ℓ-1)
                    if n == 2 && m ==1 && ℓ >= 2
                        M₂_index_α₁₀ = M₂_index_α₁₀ +1 
                    end
                    if m == 2 && n ==1 && ℓ >= 2
                        M₂_index_α₀₁ = M₂_index_α₀₁ +1 
                    end
                end
            end
        end
    end

    dN_temp = Vector{Any}(undef, 29) 
    dM₁ = convert( Int64 ,  maximum(M₁_index)  )   + 2 
    dM₂ = convert( Int64 ,  maximum(M₂_index)  )   + 2 

    dN_temp[1] =  3
    dN_temp[2] =  3
    dN_temp[3] =  dM₁
    dN_temp[4] =  3
    dN_temp[5] =  dM₂
    dN_temp[6] =  2
    dN_temp[7] =  3
    dN_temp[8] =  max( dM₁ , 3 )
    dN_temp[9] =  3
    dN_temp[10] = max( dM₂ , 3 )
    dN_temp[11] = 3
    dN_temp[12] = max( dM₁ , 3 )
    dN_temp[13] = 3
    dN_temp[14] = max( dM₂ , 3 )
    dN_temp[15] = sum(order_p[1]) + 2 
    dN_temp[16] = sum(order_p[2]) + 2 
    dN_temp[17] = sum(order_p[3]) + 2 
    dN_temp[18] = sum(order_p[4]) + 2 
    dN_temp[19] =  order_w[1][1] +2
    dN_temp[20] =  max( dM₁ , order_w[2][1] +2 ) #order(w)[2][1] +2
    dN_temp[21] =  order_w[3][1] +2 
    dN_temp[22] =  max( dM₂ , order_w[4][1] +2 ) #order(w)[4][1] +2 
    dN_temp[23] =  3
    dN_temp[24] =  3
    dN_temp[25] =  max( maximum([ M₁_index_α₁₀ ; M₁_index_α₀₁ ; M₂_index_α₁₀ ; M₂_index_α₀₁ ]) + 2 , 3 )
    dN_temp[26] =  3
    dN_temp[27] =  max( dM₁ -1 , 3 )
    dN_temp[28] =  3
    dN_temp[29] =  max( dM₂ -1 , 3 )

    dN = N*dN_temp
    return dN
end

function F_fft!(F_fft::Vector, space_X_int, space_F_int , dN_ftt, X_cheb, η_cheb, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂)
    iter = ProgressBar(1:29)
    set_description(iter, "    - F:")
    for i in iter
        F_fft[i] = map(A -> interval.(A), [zeros(ComplexF64, space_F_int[i] ) for _ ∈ 1:dN_ftt[i]])
        X_temp = cheb2grid(X_cheb, dN_ftt[i])
        η_temp = cheb2grid(η_cheb, dN_ftt[i])
        for j = 1:dN_ftt[i]
            F_fft[i][j] = Fᵢ!(F_fft[i][j], Sequence( space_X_int, X_temp[j]), κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Sequence( space_X_int, X_temp[j]) , Sequence( space_X_int, η_temp[j]) , i )
        end
    end
    return F_fft
end


function F_grid!(F_grid, space_X_int, space_F_int , dN, X_Taylor, η_Taylor , κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂)
    for i = 1:29
        F_grid[i] = map(A -> interval.(A), [zeros(ComplexF64, space_F_int[i] ) for _ ∈ 1:dN[i]+1])
        X_temp = taylor2grid_int(X_Taylor, dN[i])
        η_temp = taylor2grid_int(η_Taylor, dN[i])
        for j = 1:dN[i]+1
            F_grid[i][j] = Fᵢ!(F_grid[i][j], Sequence( space_X_int, X_temp[j][:]), κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Sequence( space_X_int, X_temp[j][:]) , Sequence( space_X_int, η_temp[j][:]) , i )
        end
    end
    return F_grid
end

function AF_fft!(AF_fft::Vector{<:Vector}, A_cheb::Vector, F_fft::Vector, space_X_int::CartesianProduct, space_F_int::CartesianProduct, dN_ftt::Vector{Int64})
    iter = ProgressBar(1:29)
    set_description(iter, "    - AF:")
    for i in iter
        test = cheb2grid(A_cheb[i], dN_ftt[i])
        AF_fft_temp = map(A -> LinearOperator( space_F_int[i], space_X_int,A), [ test[j]  for j ∈ 1:dN_ftt[i]]).*F_fft[i] 
        AF_fft[i] = map(A -> A[:], [AF_fft_temp[j] for j ∈ 1:dN_ftt[i]])
    end   
    return AF_fft 
end
function AF_grid!(AF_grid, A_Taylor, F_grid, space_X_int, space_F_int, dN)
    for i = 1:29
        test = taylor2grid_int_A_col(A_Taylor[i], dN[i] )
        AF_grid_temp = map(A -> LinearOperator( space_F_int[i], space_X_int,A), [ test[j]  for j ∈ 1:dN[i]+1]).*F_grid[i] 
        AF_grid[i] = map(A -> A[:], [AF_grid_temp[j] for j ∈ 1:dN[i]+1])
    end   
    return AF_grid 
end
function AF_cheb!( AF_fft, dN)
    AFᵢ_cheb = Vector{Any}( undef, 29 )
    for i = 1:29
        AFᵢ_cheb[i] =  grid2cheb(AF_fft[i], dN[i])
    end
    AF_cheb = copy(AFᵢ_cheb[1])
    for i = 2:29
        AF_cheb = AF_cheb + AFᵢ_cheb[i] 
    end
    return norm.(AF_cheb, 1)
end
function AF_cheb2!( AF_fft, dN)
    AFᵢ_cheb = Vector{Any}( undef, 29 )
    for i = 1:29
        AFᵢ_cheb[i] =  grid2cheb(AF_fft[i], dN[i])
    end
    AF_cheb = copy(AFᵢ_cheb[1])
    for i = 2:29
        AF_cheb = AF_cheb + AFᵢ_cheb[i] 
    end
    return norm.(AF_cheb, 1)
end
function F_ext_fft!( F_ext_fft, X_cheb, η_cheb, space_X_int, space_F_int_ext , dN_ftt, κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂)
    iter = ProgressBar(1:29)
    set_description(iter, "    - ")
    for i in iter 
        F_ext_fft[i] = map(A -> interval.(A), [zeros(ComplexF64, space_F_int_ext[i] ) for _ ∈ 1:dN_ftt[i]])
        X_temp = cheb2grid(X_cheb, dN_ftt[i])
        η_temp = cheb2grid(η_cheb, dN_ftt[i])
        for j = 1:dN_ftt[i]
            F_ext_fft[i][j] = Fᵢ!(F_ext_fft[i][j], Sequence( space_X_int, X_temp[j]), κ, ξ₁ᵏ ,ξ₂ᵏ, star, lengthᵥ, M₁, M₂, Sequence( space_X_int, X_temp[j]) , Sequence( space_X_int, η_temp[j]) , i )
        end
    end
    return F_ext_fft
end

function Y₀_tails_fft!( Y₀_tails_fft, λ_cheb,λ₁_cheb,λ₂_cheb, dN, dN_ftt,space_F ,w ,a,p ,space_F_int_ext , Λ_tails) 
    for i = 1:29
       Λ_tails[i] = Vector{Any}(undef, dN_ftt[i])
        λ_temp = cheb2grid(λ_cheb, dN_ftt[i])
        λ₁_temp = cheb2grid(λ₁_cheb, dN_ftt[i])
        λ₂_temp = cheb2grid(λ₂_cheb, dN_ftt[i])
        for j = 1:dN_ftt[i]
            if i in [1 ; 6 ; 15  ; 16 ; 17 ; 18 ; 23 ; 24 ; 25]
                Λ_tails[i][j] = ExactReal(0)
                Y₀_tails_fft[i][j] .= ExactReal(0)
            elseif i in [2 ; 4]
                Λ_tails[i][j] = interval(1) ./ (interval.(order(space_F[i])) .+ interval(1))
                Y₀_tails_fft[i][j][1:order(space_F[i])] .= ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            elseif i in [ 3 ; 5]
                Λ_tails[i][j] = interval(1) ./ (interval.(order(space_F[i])) .+ interval(1))
                Y₀_tails_fft[i][j][0:order(space_F[i])] .= ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            elseif i in [7 ; 8 ; 9 ; 10]
                Λ_tails[i][j] =  norm(interval(1)/( λ_temp[j][1] - interval(1im*(order(space_F[i])+1)) )) 
                Y₀_tails_fft[i][j][-order(space_F[i]):order(space_F[i])] .= ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            elseif i in [11 ; 12 ; 13 ; 14]
                n = ExactReal.(range(-order(w)[i-10][1]-1,order(w)[i-10][1]+1))
                m = ExactReal.(range(0,order(w)[i-10][2]+1))
                C₁ = maximum(norm.(interval(1) ./ (interval(1im)*(interval(order(w)[i-10][1]+1)) .+ λ_temp[j][1]*m)))
                C₂ = maximum(norm.(interval(1) ./ (-interval(1im)*(interval(order(w)[i-10][1]+1)) .+ λ_temp[j][1]*m)))
                Λ_tails[i][j] = maximum( [ C₁ , C₂ ]  ) 
                Y₀_tails_fft[i][j][ (-order(space_F[i])[1]:order(space_F[i])[1], 0:order(space_F[i])[2]) ] .=ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            elseif i in [19; 20; 21; 22]
                Λ_tails[i][j] = (interval(1) ./ ( interval(2)*interval((order(a)[i-18]+1))))
                Y₀_tails_fft[i][j][0:order(space_F[i])] .= ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            elseif i in [26; 27; 28; 29]
                n = ExactReal.(range(0,order(p)[i-25][1]+1))
                m = ExactReal.(range(0,order(p)[i-25][2]+1))
                C₁ = maximum(norm.(interval(1) ./ (λ₁_temp[j][1]*n .+ λ₂_temp[j][1]*interval(order(p)[i-25][2]+1)) ))
                C₂ = maximum(norm.(interval(1) ./ (λ₁_temp[j][1]*interval(order(p)[i-25][1]+1) .+ λ₂_temp[j][1]*m )))
                Λ_tails[i][j] = maximum( [  C₁,C₂]  )
                Y₀_tails_fft[i][j][ (0:order(space_F[i])[1], 0:order(space_F[i])[2]) ] .= ExactReal(0)
                Y₀_tails_fft[i][j] =  Λ_tails[i][j].* Y₀_tails_fft[i][j]
            end
        end
    end
    Y₀_tails_cheb =  zeros( Complex{Interval{Float64}}, space_F_int_ext  )
    for i = 1:29
        Y₀_tails_temp = Vector{Vector}(undef, dN_ftt[i] )
        for j = 1:dN_ftt[i]
            Y₀_tails_temp[j] =  Y₀_tails_fft[i][j][:]
        end
        component( Y₀_tails_cheb ,i)[:] =  norm.(  grid2cheb(Y₀_tails_temp, dN[i]) , 1  ) 
    end

    return Y₀_tails_cheb 
end

function DN!( N , M₁ , M₂ , order_p , order_w, κ)
    M₁_index =  M₁.*0
    M₁_index_α₁₀ =  0
    M₁_index_α₀₁ =  0
    for n in axes(M₁,1)
        for m in axes(M₁,2)
            for ℓ in axes(M₁,3)
                if M₁[n,m,ℓ] != 0
                    M₁_index[n,m,ℓ] = (n-1) + (m-1) + (ℓ-1)
                    if n == 2 && m ==1 && ℓ >= 2
                        M₁_index_α₁₀ = M₁_index_α₁₀ +1 
                    end
                    if m == 2 && n ==1 && ℓ >= 2
                        M₁_index_α₀₁ = M₁_index_α₀₁ +1 
                    end
                end
            end
        end
    end    
    M₂_index =  M₂.*0
    M₂_index_α₁₀ =  0
    M₂_index_α₀₁ =  0
    for n in axes(M₂,1)
        for m in axes(M₂,2)
            for ℓ in axes(M₂,3)
                if M₂[n,m,ℓ] != 0
                    M₂_index[n,m,ℓ] = (n-1) + (m-1) + (ℓ-1)
                    if n == 2 && m ==1 && ℓ >= 2
                        M₂_index_α₁₀ = M₂_index_α₁₀ +1 
                    end
                    if m == 2 && n ==1 && ℓ >= 2
                        M₂_index_α₀₁ = M₂_index_α₀₁ +1 
                    end
                end
            end
        end
    end

    dM₁ = convert( Int64 ,  maximum(M₁_index)  )   + 2 
    dM₂ = convert( Int64 ,  maximum(M₂_index)  )   + 2 

    DN_temp = zeros(Int64,29,29)
    # i = 1
    DN_temp[1,:] .= 2 
    # i = 2
    DN_temp[2,2] = 1 
    DN_temp[2,3] = 2 
    DN_temp[2,24] = 2 
    # i = 3
    DN_temp[3,1] = dM₁-1
    DN_temp[3,2] = dM₁-1
    DN_temp[3,3] = 1 
    DN_temp[3,4] = dM₁-1
    DN_temp[3,24] = dM₁-1
    # i = 4
    DN_temp[4,4] = 1 
    DN_temp[4,5] = 2 
    DN_temp[4,24] = 2 
    # i = 5
    DN_temp[5,1] = dM₂-1
    DN_temp[5,2] = dM₂-1
    DN_temp[5,4] = dM₂-1
    DN_temp[5,5] = 1 
    DN_temp[5,24] = dM₂-1
    # i = 6
    DN_temp[6,6+κ] = 1
    # i = 7
    DN_temp[ 7 , 6 ] = 2
    DN_temp[ 7 , 7 ] = 2
    DN_temp[ 7 , 8 ] = 2
    DN_temp[ 7 , 24 ] = 2
    # i = 8
    DN_temp[ 8 , 1 ] = dM₁ - 1
    DN_temp[ 8 , 2 ] = dM₁ - 1
    DN_temp[ 8 , 4 ] = dM₁ - 1
    DN_temp[ 8 , 6 ] = 2
    DN_temp[ 8 , 7 ] = dM₁ - 1
    DN_temp[ 8 , 8 ] = 2
    DN_temp[ 8 , 9 ] = dM₁ - 1
    DN_temp[ 8 , 24 ] = dM₁ - 1
    # i = 9
    DN_temp[ 9 , 6 ] = 2
    DN_temp[ 9 , 9 ] = 2
    DN_temp[ 9 , 10 ] = 2
    DN_temp[ 9 , 24 ] = 2
    # i = 10
    DN_temp[ 10 , 1 ] = dM₂ - 1
    DN_temp[ 10 , 2 ] = dM₂ - 1
    DN_temp[ 10 , 4 ] = dM₂ - 1
    DN_temp[ 10 , 6 ] = 2
    DN_temp[ 10 , 7 ] = dM₂ - 1
    DN_temp[ 10 , 9 ] = 2
    DN_temp[ 10 , 10 ] = dM₂ - 1
    DN_temp[ 10 , 24 ] = dM₂ - 1
    # i = 11
    DN_temp[ 11 , 2  ] = 1
    DN_temp[ 11 , 6  ] = 2
    DN_temp[ 11 , 7  ] = 1
    DN_temp[ 11 , 11 ] = 2
    DN_temp[ 11 , 12 ] = 2
    DN_temp[ 11 , 24 ] = 2
    # i = 12
    DN_temp[ 12 , 1  ] = dM₁-1
    DN_temp[ 12 , 3  ] = 1
    DN_temp[ 12 , 6  ] = 2
    DN_temp[ 12 , 8  ] = 1
    DN_temp[ 12 , 12 ] = 2
    DN_temp[ 12 , 11 ] = dM₁-1
    DN_temp[ 12 , 13 ] = dM₁-1
    DN_temp[ 12 , 24 ] = dM₁-1
    # i = 13
    DN_temp[ 13 , 4  ] = 1
    DN_temp[ 13 , 6  ] = 2
    DN_temp[ 13 , 9  ] = 1
    DN_temp[ 13 , 13 ] = 2
    DN_temp[ 13 , 14 ] = 2
    DN_temp[ 13 , 24 ] = 2
    # i = 14
    DN_temp[ 14 , 1  ] = dM₂-1
    DN_temp[ 14 , 5  ] = 1
    DN_temp[ 14 , 6  ] = 2
    DN_temp[ 14 , 10  ] = 1
    DN_temp[ 14 , 14 ] = 2
    DN_temp[ 14 , 11 ] = dM₂-1
    DN_temp[ 14 , 13 ] = dM₂-1
    DN_temp[ 14 , 24 ] = dM₂-1
    # i = 15
    DN_temp[15 , 18] = sum(order_p[1]) + 1
    DN_temp[15 , 19] = sum(order_p[1]) + 1
    DN_temp[15 , 20] = 1
    DN_temp[15 , 26] = sum(order_p[1]) + 1
    # i = 16
    DN_temp[16 , 18] = sum(order_p[1]) + 1
    DN_temp[16 , 19] = sum(order_p[1]) + 1
    DN_temp[16 , 21] = 1
    DN_temp[16 , 27] = sum(order_p[1]) + 1
    # i = 17
    DN_temp[17 , 18] = sum(order_p[1]) + 1
    DN_temp[17 , 19] = sum(order_p[1]) + 1
    DN_temp[17 , 22] = 1
    DN_temp[17 , 28] = sum(order_p[1]) + 1
    # i = 18
    DN_temp[18 , 18] = sum(order_p[1]) + 1
    DN_temp[18 , 19] = sum(order_p[1]) + 1
    DN_temp[18 , 23] = 1
    DN_temp[18 , 29] = sum(order_p[1]) + 1
    # i = 19
    DN_temp[ 19 , 11 ] = order_w[1][1] +1
    DN_temp[ 19 , 15 ] = 2
    DN_temp[ 19 , 16 ] = order_w[1][1] +1
    DN_temp[ 19 , 17 ] = order_w[1][1] +1
    DN_temp[ 19 , 20 ] = 1
    DN_temp[ 19 , 21 ] = 2
    # i = 20
    DN_temp[ 20 , 12 ] = order_w[2][1] +1
    DN_temp[ 20 , 15 ] = 2
    DN_temp[ 20 , 16 ] = order_w[2][1] +1
    DN_temp[ 20 , 17 ] = order_w[2][1] +1
    DN_temp[ 20 , 21 ] = 1
    DN_temp[ 20 , 20 ] = dM₁-1
    DN_temp[ 20 , 22 ] = dM₁-1
    DN_temp[ 20 ,  1 ] = dM₁-1
    # i = 21
    DN_temp[ 21 , 13 ] = order_w[3][1] +1
    DN_temp[ 21 , 15 ] = 2
    DN_temp[ 21 , 16 ] = order_w[3][1] +1
    DN_temp[ 21 , 17 ] = order_w[3][1] +1
    DN_temp[ 21 , 22 ] = 1
    DN_temp[ 21 , 23 ] = 2
    # i = 22
    DN_temp[ 22 , 14 ] = order_w[4][1] +1
    DN_temp[ 22 , 15 ] = 2
    DN_temp[ 22 , 16 ] = order_w[4][1] +1
    DN_temp[ 22 , 17 ] = order_w[4][1] +1
    DN_temp[ 22 , 23 ] = 1
    DN_temp[ 22 , 20 ] = dM₂-1
    DN_temp[ 22 , 22 ] = dM₂-1
    DN_temp[ 22 ,  1 ] = dM₂-1
    # i = 23
    DN_temp[ 23 , 18 ] = 2
    DN_temp[ 23 , 19 ] = 2
    # i = 24
    DN_temp[ 24 , 16 ] = 2
    DN_temp[ 24 , 17 ] = 2
    # i = 25
    DN_temp[ 25 , 1 ] = max( maximum([ M₁_index_α₁₀ ; M₁_index_α₀₁ ; M₂_index_α₁₀ ; M₂_index_α₀₁ ]) + 2 , 3 ) -1
    DN_temp[ 25 , 25 ] = max( maximum([ M₁_index_α₁₀ ; M₁_index_α₀₁ ; M₂_index_α₁₀ ; M₂_index_α₀₁ ]) + 2 , 3 ) -1
    # i = 26
    DN_temp[ 26 , 25 ] = 2
    DN_temp[ 26 , 26 ] = 2
    DN_temp[ 26 , 27 ] = 1
    # i = 27
    DN_temp[ 27 , 25 ] = 2
    DN_temp[ 27 , 27 ] = 2
    DN_temp[ 27 , 26 ] = dM₁ -1
    DN_temp[ 27 , 28 ] = dM₁ -1
    DN_temp[ 27 , 1  ] = dM₁ -1
    # i = 28
    DN_temp[ 28 , 25 ] = 2
    DN_temp[ 28 , 28 ] = 2
    DN_temp[ 28 , 29 ] = 1
    # i = 29
    DN_temp[ 29 , 25 ] = 2
    DN_temp[ 29 , 29 ] = 2
    DN_temp[ 29 , 26 ] = dM₁ -1
    DN_temp[ 29 , 28 ] = dM₁ -1
    DN_temp[ 29 , 1  ] = dM₁ -1
    DN = N*(DN_temp)
    return DN
end

function Aᵢⱼ_cheb!(Aᵢⱼ_cheb, space_X , space_F, N, N_fft , X_fft, κ, star, lengthᵥ, M₁, M₂ , η_fft)
    Aᵢⱼ_fft = Matrix{Vector}( undef, 29 , 29 )
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N_fft]
    for ℓ  = 1:N_fft
        DF_ℓ = inv(DF_all!(DF[ℓ], X_fft[ℓ], κ, star, lengthᵥ, M₁, M₂ , η_fft[ℓ] ))
        for i = 1:29
            for j = 1:29
                if ℓ == 1
                    Aᵢⱼ_fft[i,j] = Vector{Matrix}(undef, N_fft)
                end
                Aᵢⱼ_fft[i,j][ℓ] = component( DF_ℓ , i , j)[:,:]
            end
        end
    end
    for i = 1:29
        for j = 1:29
            Aᵢⱼ_cheb[i,j] = map(A -> interval.(A), grid2cheb(Aᵢⱼ_fft[i,j], N))
        end
    end
    return Aᵢⱼ_cheb
end

function Aᵢⱼ_cheb_v2!(Aᵢⱼ_cheb, space_X , space_F, N, N_fft , X_fft, κ, star, lengthᵥ, M₁, M₂ , η_fft)
    Aᵢⱼ_fft = Matrix{Vector}( undef, 29 , 29 )
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N_fft]
    for ℓ  = 1:N_fft
        DF_ℓ = inv(DF_all!(DF[ℓ], X_fft[ℓ], κ, star, lengthᵥ, M₁, M₂ , η_fft[ℓ] ))
        for i = 1:29
            for j = 1:29
                if ℓ == 1
                    Aᵢⱼ_fft[i,j] = Vector{Matrix}(undef, N_fft)
                end
                Aᵢⱼ_fft[i,j][ℓ] = component( DF_ℓ , i , j)[:,:]
            end
        end
    end
    for i = 1:29
        for j = 1:29
            Aᵢⱼ_cheb[i,j] =  grid2cheb(Aᵢⱼ_fft[i,j], N)
        end
    end
    return Aᵢⱼ_cheb
end

function Aᵢⱼ_fft_v2!(Aᵢⱼ_fft, space_X , space_F, N, N_fft , X_fft, κ, star, lengthᵥ, M₁, M₂ , η_fft)
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N_fft]
    for ℓ  = 1:N_fft
        Aᵢⱼ_fft[ℓ]= inv(DF_all!(DF[ℓ], X_fft[ℓ], κ, star, lengthᵥ, M₁, M₂ , η_fft[ℓ] ))[:,:]
    end
    return Aᵢⱼ_fft
end


function Aᵢⱼ_cheb_v2_int!(Aᵢⱼ_cheb, space_X , space_F, N, N_fft , X_grid, κ, star, lengthᵥ, M₁, M₂ , η_grid)
    Aᵢⱼ_fft = Matrix{Vector}( undef, 29 , 29 )
    DF = [zeros(ComplexF64, space_X , space_F) for _ ∈ 1:N+1]
    for ℓ  = 1:N+1
        DF_ℓ = inv(DF_all!(DF[ℓ], X_grid[ℓ], κ, star, lengthᵥ, M₁, M₂ , η_grid[ℓ] ))
        for i = 1:29
            for j = 1:29
                if ℓ == 1
                    Aᵢⱼ_fft[i,j] = Vector{Matrix}(undef, N+1)
                end
                Aᵢⱼ_fft[i,j][ℓ] = component( DF_ℓ , i , j)[:,:]
            end
        end
    end
    for i = 1:29
        for j = 1:29
            Aᵢⱼ_cheb[i,j] =  dft_grid2cheb(Aᵢⱼ_fft[i,j])
        end
    end
    return Aᵢⱼ_cheb
end


function Z₁_mat!(Z₁_, X_cheb, DN, DN_ftt, Λ , space_X_int, M₁_int, M₂_int , 𝑿ᵧ , 𝑿ᵥ , 𝑿_w , 𝑿ₚ , 𝑿ₐ  ,νₚ , νₐ ,ν_w , Aᵢⱼ_cheb,space_X_int_ext,space_F_int, νᵥ ,νᵧ)
    Z₁_mat = Matrix{Vector}(undef, 29 , 29)
    X_cheb_sequence = Sequence( space_X_int, X_cheb)
    σ₁_cheb = component(X_cheb_sequence , 18)[1]
    σ₂_cheb = component(X_cheb_sequence , 19)[1]
    norm_z = sqrt( norm(σ₁_cheb^2 + σ₂_cheb^2,ExactReal(1)) )

    θ₁_cheb = component(X_cheb_sequence , 16)[1]
    θ₂_cheb = component(X_cheb_sequence , 17)[1]
    norm_y = sqrt( norm(θ₁_cheb^2 + θ₂_cheb^2,ExactReal(1)) )

    norm_T = ExactReal(1)/νₐ + ExactReal(2)*νₐ
    
    order_w₁ = order(space_X_int_ext[11])
    order_w₂ = order(space_X_int_ext[12])
    order_w₃ = order(space_X_int_ext[13])
    order_w₄ = order(space_X_int_ext[14])

    order_a₁ = order(space_X_int_ext[20])
    order_a₂ = order(space_X_int_ext[21])
    order_a₃ = order(space_X_int_ext[22])
    order_a₄ = order(space_X_int_ext[23])

    order_p₁ = order(space_X_int_ext[26])
    order_p₂ = order(space_X_int_ext[27])
    order_p₃ = order(space_X_int_ext[28])
    order_p₄ = order(space_X_int_ext[29])

    iter = ProgressBar(1:29)
    set_description(iter, "    - ")
    for i in iter
        for j = 1:29
            Z₁_mat[i,j] = zeros( eltype(Z₁_) , DN_ftt[i,j] )
            if DN_ftt[i,j] > 1
                X_temp = cheb2grid( X_cheb   , DN_ftt[i,j]  )
                local Λ_tails = Λ[i]
                for k = 1:DN_ftt[i,j]
                    local X_ = Sequence( space_X_int , X_temp[k]  )
                    local α_int = real(X_[1])
                    local x = x_remove_complex!(component(X_,2:29))
                    local γ,λ,v,w,L,θ₀,σ,a,ω_int,eigenpairs,p = x2var!(x)
                    local γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
                    local γ_ext = γ_extened_orientable!(γ)
                    local γ₁_ext,γ₂_ext,γ₃_ext,γ₄_ext = eachcomponent(γ_ext)
                    local v₁,v₂,v₃,v₄ = eachcomponent(v)
                    local w₁, w₂, w₃, w₄ = eachcomponent(w)
                    local p₁, p₂, p₃, p₄ = eachcomponent(p)
                    local a₁, a₂, a₃, a₄ = eachcomponent(a)

                    if i == 2 && j == 3 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int) end # ✓
                    if i == 4 && j == 5 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int) end # ✓

                    if i == 7 && j == 8 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int) end # ✓
                    if i == 9 && j == 10 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int) end # ✓

                    if i == 11 && j == 2 Z₁_mat[i,j][k] = Λ_tails end # Missing
                    if i == 13 && j == 4 Z₁_mat[i,j][k] = Λ_tails end # Missing

                    if i == 11 && j == 7 Z₁_mat[i,j][k] = ν_w*Λ_tails end # Missing
                    if i == 13 && j == 9 Z₁_mat[i,j][k] = ν_w*Λ_tails end # Missing

                    if i == 11 && j == 12 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int)  end # Missing
                    if i == 13 && j == 14 Z₁_mat[i,j][k] = Λ_tails*abs(ω_int) end # Missing

                    if i == 12 && j == 3 Z₁_mat[i,j][k] = Λ_tails end # Missing
                    if i == 14 && j == 5 Z₁_mat[i,j][k] = Λ_tails end # Missing

                    if i == 12 && j == 8 Z₁_mat[i,j][k] = ν_w*Λ_tails end # Missing
                    if i == 14 && j == 10 Z₁_mat[i,j][k] = ν_w*Λ_tails end # Missing

                    if i == 19 && j == 21 Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T end
                    if i == 21 && j == 23 Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T end

                    if i == 3 && j == 2  Z₁_Γ₂_ᵧ₁ = zeros(eltype(γ₁), space(γ₁) ) end
                    if i == 3 && j == 4  Z₁_Γ₂_ᵧ₃ = zeros(eltype(γ₃), space(γ₃) ) end
                    if i == 3 && j == 1  Z₁_Γ₂_α = zeros(eltype(γ₃), space(γ₃) ) end
                    if i == 3 && j == 24  Z₁_Γ₂_ω = zeros(eltype(γ₃), space(γ₃) ) end

                    if i == 8 && j == 7  Z₁_V₂_ᵥ₁ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 8 && j == 9  Z₁_V₂_ᵥ₃ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 8 && j == 2  Z₁_V₂_ᵧ₁ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 8 && j == 4  Z₁_V₂_ᵧ₃ = zeros(eltype(γ₃), space(v₁) ) end
                    if i == 8 && j == 1  Z₁_V₂_α = zeros(eltype(v₁), space(v₁) ) end
                    if i == 8 && j == 24  Z₁_V₂_ω = zeros(eltype(v₁), space(v₁) ) end
                    
                    if i == 12 && j == 11  Z₁_W₂_w₁ = zeros(eltype(w₁), space(w₁) ) end
                    if i == 12 && j == 13 Z₁_W₂_w₃ = zeros(eltype(w₃), space(w₃) ) end
                    if i == 12 && j == 1  Z₁_W₂_α = zeros(eltype(w₃), space(w₃) ) end
                    if i == 12 && j == 24  Z₁_W₂_ω = zeros(eltype(w₃), space(w₃) ) end

                    if i == 27 && j == 26  Z₁_P₂_ₚ₁ = zeros(eltype(p₁), space(p₁) ) end
                    if i == 27 && j == 28  Z₁_P₂_ₚ₃ = zeros(eltype(p₃), space(p₃) ) end
                    if i == 27 && j == 1  Z₁_P₂_α = zeros(eltype(p₃), space(p₃) ) end

                    if i == 20 && j == 20  Z₁_G₂_ₐ₁ = zeros(eltype(a₁), space(a₁) ) end
                    if i == 20 && j == 22  Z₁_G₂_ₐ₃ = zeros(eltype(a₃), space(a₃) ) end
                    if i == 20 && j == 1  Z₁_G₂_α = zeros(eltype(a₃), space(a₃) ) end
                    if i == 20 && j == 15  Z₁_G₂_L = zeros(eltype(a₃), space(a₃) ) end

                    for n in axes(M₁_int,1)
                        for m in axes(M₁_int,2)
                            for ℓ in axes(M₁_int,3)
                                if M₁_int[n,m,ℓ] != 0 
                                    if i == 20 && j == 15  Z₁_G₂_L = Z₁_G₂_L - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * a₁^(n-1)*a₃^(m-1) end# ✓
                                    if i == 3 && j == 24  Z₁_Γ₂_ω = Z₁_Γ₂_ω - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)) end # ✓
                                    if i == 12 && j == 24  Z₁_W₂_ω = Z₁_W₂_ω - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*(w₁^(n-1)*w₃^(m-1)) end# ✓
                                    if ℓ >= 2
                                        if i == 3 && j == 1  Z₁_Γ₂_α = Z₁_Γ₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₁_int[n,m,ℓ] * γ₁^(n-1)*γ₃^(m-1) end # ✓
                                        if i == 12 && j == 1  Z₁_W₂_α = Z₁_W₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₁_int[n,m,ℓ] * w₁^(n-1)*w₃^(m-1) end# ✓
                                        if i == 27 && j == 1  Z₁_P₂_α = Z₁_P₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*M₁_int[n,m,ℓ] * p₁^(n-1)*p₃^(m-1) end
                                        if i == 20 && j == 1  Z₁_G₂_α = Z₁_G₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*M₁_int[n,m,ℓ] * a₁^(n-1)*a₃^(m-1) end# ✓
                                    end
                                    if n >= 2
                                        if i == 3 && j == 2  Z₁_Γ₂_ᵧ₁ = Z₁_Γ₂_ᵧ₁ - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(n-1)*(γ₁^(n-2)*γ₃^(m-1)) end# ✓
                                        if i == 12 && j == 11  Z₁_W₂_w₁ = Z₁_W₂_w₁ - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(n-1)*(w₁^(n-2)*w₃^(m-1)) end# ✓
                                        if i == 27 && j == 26  Z₁_P₂_ₚ₁ = Z₁_P₂_ₚ₁ - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*ExactReal(n-1)*p₁^(n-2)*p₃^(m-1) end
                                        if i == 20 && j == 20  Z₁_G₂_ₐ₁ = Z₁_G₂_ₐ₁ - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*ExactReal(n-1)*a₁^(n-2)*a₃^(m-1) end# ✓
                                    end
                                    if m >= 2
                                        if i == 3 && j == 4  Z₁_Γ₂_ᵧ₃ = Z₁_Γ₂_ᵧ₃ - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(m-1)*(γ₁^(n-1)*γ₃^(m-2)) end # ✓
                                        if i == 12 && j == 13 Z₁_W₂_w₃ = Z₁_W₂_w₃ - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(m-1)*(w₁^(n-1)*w₃^(m-2)) end # ✓
                                        if i == 27 && j == 28  Z₁_P₂_ₚ₃ = Z₁_P₂_ₚ₃ - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*ExactReal(m-1)*p₁^(n-1)*p₃^(m-2) end
                                        if  i == 20 && j == 22  Z₁_G₂_ₐ₃ = Z₁_G₂_ₐ₃ - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ]*ExactReal(m-1)*a₁^(n-1)*a₃^(m-2) end# ✓
                                    end
                                    if n-2 >= 0
                                        if i == 8 && j == 7  Z₁_V₂_ᵥ₁ = Z₁_V₂_ᵥ₁  - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)  end# ✓
                                    end
                                    if m-2 >= 0
                                        if i == 8 && j == 9  Z₁_V₂_ᵥ₃ = Z₁_V₂_ᵥ₃  - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)  end# ✓
                                    end
                                    if n-2 >= 0 && m-2>= 0
                                        if i == 8 && j == 24  Z₁_V₂_ω = Z₁_V₂_ω - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * (ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ) end# ✓
                                    elseif n-2 >= 0
                                        if i == 8 && j == 24  Z₁_V₂_ω = Z₁_V₂_ω - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ end# ✓
                                    elseif m-2 >= 0
                                        if i == 8 && j == 24  Z₁_V₂_ω = Z₁_V₂_ω - α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ end# ✓
                                    end
                                    if n-3 >= 0 
                                        if i == 8 && j == 2  Z₁_V₂_ᵧ₁ = Z₁_V₂_ᵧ₁  -ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(n-1)*ExactReal(n-2)*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ end# ✓
                                    end
                                    if n-2 >= 0 && m-2 >= 0
                                        if i == 8 && j == 2  Z₁_V₂_ᵧ₁ = Z₁_V₂_ᵧ₁  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ end# ✓
                                        if i == 8 && j == 4  Z₁_V₂_ᵧ₃ = Z₁_V₂_ᵧ₃  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ end# ✓
                                    end
                                    if m-3 >=0
                                        if i == 8 && j == 4  Z₁_V₂_ᵧ₃ = Z₁_V₂_ᵧ₃  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(m-2)*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ end# ✓
                                    end
                                    if n-2 >= 0 && m-2>= 0 && ℓ >= 2
                                        if i == 8 && j == 1  Z₁_V₂_α = Z₁_V₂_α- ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₁_int[n,m,ℓ] * (ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ) end# ✓
                                    elseif n-2 >= 0 && ℓ >= 2
                                        if i == 8 && j == 1  Z₁_V₂_α = Z₁_V₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₁_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ end# ✓
                                    elseif m-2 >= 0 && ℓ >= 2
                                        if i == 8 && j == 1  Z₁_V₂_α = Z₁_V₂_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₁_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ end# ✓
                                    end
                                end
                            end
                        end
                    end
                    
                    if i == 3 && j == 2  Z₁_mat[i,j][k]  =  Λ_tails*norm(Z₁_Γ₂_ᵧ₁, 𝑿ᵧ )  end
                    if i == 3 && j == 4  Z₁_mat[i,j][k]   =  Λ_tails*norm(Z₁_Γ₂_ᵧ₃, 𝑿ᵧ ) end
                    if i == 3 && j == 24  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_Γ₂_ω, 𝑿ᵧ ) end
                    if i == 3 && j == 1  Z₁_mat[i,j][k]  =  Λ_tails*norm(Z₁_Γ₂_α, 𝑿ᵧ ) end

                    if i == 8 && j == 2  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_ᵧ₁, 𝑿ᵥ ) end
                    if i == 8 && j == 4  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_ᵧ₃, 𝑿ᵥ ) end
                    if i == 8 && j == 24  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_ω, 𝑿ᵥ )  end
                    if i == 8 && j == 1  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_α, 𝑿ᵥ ) end
                    if i == 8 && j == 7  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_ᵥ₁, 𝑿ᵥ ) end
                    if i == 8 && j == 9  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₂_ᵥ₃, 𝑿ᵥ ) end

                    if i == 12 && j == 11  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₂_w₁, 𝑿_w ) end
                    if i == 12 && j == 13 Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₂_w₃, 𝑿_w ) end
                    if i == 12 && j == 24  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₂_ω, 𝑿_w )  end
                    if i == 12 && j == 1  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₂_α, 𝑿_w ) end

                    if i == 27 && j == 26  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₂_ₚ₁, 𝑿ₚ ) end
                    if i == 27 && j == 28  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₂_ₚ₃, 𝑿ₚ ) end
                    if i == 27 && j == 1  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₂_α, 𝑿ₚ ) end

                    if i == 20 && j == 20  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₂_ₐ₁, 𝑿ₐ )  end
                    if i == 20 && j == 22  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₂_ₐ₃, 𝑿ₐ ) end
                    if i == 20 && j == 1  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₂_α, 𝑿ₐ ) end
                    if i == 20 && j == 15  Z₁_mat[i,j][k] = Λ_tails/ExactReal(4)*norm_T*norm(Z₁_G₂_L, 𝑿ₐ ) end

                    if i == 5 && j == 2  Z₁_Γ₄_ᵧ₁ = zeros(eltype(γ₁), space(γ₁) ) end
                    if i == 5 && j == 4  Z₁_Γ₄_ᵧ₃ = zeros(eltype(γ₃), space(γ₃) ) end
                    if i == 5 && j == 1  Z₁_Γ₄_α = zeros(eltype(γ₃), space(γ₃) ) end
                    if i == 5 && j == 24  Z₁_Γ₄_ω = zeros(eltype(γ₃), space(γ₃) ) end

                    if i == 10 && j == 7  Z₁_V₄_ᵥ₁ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 10 && j == 9  Z₁_V₄_ᵥ₃ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 10 && j == 2  Z₁_V₄_ᵧ₁ = zeros(eltype(v₁), space(v₁) ) end
                    if i == 10 && j == 4  Z₁_V₄_ᵧ₃ = zeros(eltype(γ₃), space(v₁) ) end
                    if i == 10 && j == 1   Z₁_V₄_α = zeros(eltype(v₁), space(v₁) ) end
                    if i == 10 && j == 24  Z₁_V₄_ω = zeros(eltype(v₁), space(v₁) ) end

                    if i == 14 && j == 11  Z₁_W₄_w₁ = zeros(eltype(w₁), space(w₁) ) end
                    if i == 14 && j == 13  Z₁_W₄_w₃ = zeros(eltype(w₃), space(w₃) ) end
                    if i == 14 && j == 1  Z₁_W₄_α = zeros(eltype(w₃), space(w₃) ) end
                    if i == 14 && j == 24  Z₁_W₄_ω = zeros(eltype(w₃), space(w₃) ) end

                    if i == 29 && j == 26  Z₁_P₄_ₚ₁ = zeros(eltype(p₁), space(p₁) ) end
                    if i == 29 && j == 28  Z₁_P₄_ₚ₃ = zeros(eltype(p₃), space(p₃) ) end
                    if i == 29 && j == 1  Z₁_P₄_α = zeros(eltype(p₃), space(p₃) ) end

                    if i == 22 && j == 20  Z₁_G₄_ₐ₁ = zeros(eltype(a₁), space(a₁) ) end
                    if i == 22 && j == 22  Z₁_G₄_ₐ₃ = zeros(eltype(a₃), space(a₃) ) end
                    if i == 22 && j == 1  Z₁_G₄_α = zeros(eltype(a₃), space(a₃) ) end
                    if i == 22 && j == 15  Z₁_G₄_L = zeros(eltype(a₃), space(a₃) ) end

                    for n  in axes(M₂_int,1)
                        for m in axes(M₂_int,2)
                            for ℓ in axes(M₂_int,3)
                                if M₂_int[n,m,ℓ] != 0
                                    if i == 22 && j == 15  Z₁_G₄_L = Z₁_G₄_L - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * a₁^(n-1)*a₃^(m-1)end
                                    if i == 5 && j == 24  Z₁_Γ₄_ω = Z₁_Γ₄_ω - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*(γ₁^(n-1)*γ₃^(m-1)) end # ✓
                                    if i == 14 && j == 24  Z₁_W₄_ω = Z₁_W₄_ω - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*(w₁^(n-1)*w₃^(m-1)) end
                                    if ℓ >= 2
                                        if i == 5 && j == 1  Z₁_Γ₄_α = Z₁_Γ₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₂_int[n,m,ℓ] * γ₁^(n-1)*γ₃^(m-1) end # ✓
                                        if i == 14 && j == 1  Z₁_W₄_α = Z₁_W₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₂_int[n,m,ℓ] * w₁^(n-1)*w₃^(m-1) end
                                        if i == 29 && j == 1  Z₁_P₄_α = Z₁_P₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*M₂_int[n,m,ℓ] * p₁^(n-1)*p₃^(m-1) end
                                        if i == 22 && j == 1  Z₁_G₄_α = Z₁_G₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*M₂_int[n,m,ℓ] * a₁^(n-1)*a₃^(m-1) end
                                    end
                                    if n >= 2
                                        if i == 5 && j == 2  Z₁_Γ₄_ᵧ₁ = Z₁_Γ₄_ᵧ₁ - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(n-1)*(γ₁^(n-2)*γ₃^(m-1)) end # ✓
                                        if i == 14 && j == 11  Z₁_W₄_w₁ = Z₁_W₄_w₁ - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(n-1)*(w₁^(n-2)*w₃^(m-1)) end
                                        if i == 29 && j == 26  Z₁_P₄_ₚ₁ = Z₁_P₄_ₚ₁ - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*ExactReal(n-1)*p₁^(n-2)*p₃^(m-1) end
                                        if i == 22 && j == 20  Z₁_G₄_ₐ₁ = Z₁_G₄_ₐ₁ - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*ExactReal(n-1)*a₁^(n-2)*a₃^(m-1) end
                                    end
                                    if m >= 2
                                        if i == 5 && j == 4  Z₁_Γ₄_ᵧ₃ = Z₁_Γ₄_ᵧ₃ - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(m-1)*(γ₁^(n-1)*γ₃^(m-2)) end # ✓
                                        if i == 14 && j == 13  Z₁_W₄_w₃ = Z₁_W₄_w₃ - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(m-1)*(w₁^(n-1)*w₃^(m-2)) end
                                        if i == 29 && j == 28  Z₁_P₄_ₚ₃ = Z₁_P₄_ₚ₃ - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*ExactReal(m-1)*p₁^(n-1)*p₃^(m-2) end
                                        if i == 22 && j == 22  Z₁_G₄_ₐ₃ = Z₁_G₄_ₐ₃ - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ]*ExactReal(m-1)*a₁^(n-1)*a₃^(m-2) end
                                    end
                                    if n-2 >= 0
                                        if i == 10 && j == 7  Z₁_V₄_ᵥ₁ = Z₁_V₄_ᵥ₁  - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) end
                                    end
                                    if m-2 >= 0
                                        if i == 10 && j == 9  Z₁_V₄_ᵥ₃ = Z₁_V₄_ᵥ₃  - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) end
                                    end
                                    if n-2 >= 0 && m-2>= 0
                                        if i == 10 && j == 24  Z₁_V₄_ω = Z₁_V₄_ω - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * (ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ) end
                                    elseif n-2 >= 0
                                        if i == 10 && j == 24  Z₁_V₄_ω = Z₁_V₄_ω - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ end
                                    elseif m-2 >= 0
                                        if i == 10 && j == 24  Z₁_V₄_ω = Z₁_V₄_ω - α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ end
                                    end
                                    if n-3 >= 0 
                                        if i == 10 && j == 2  Z₁_V₄_ᵧ₁ = Z₁_V₄_ᵧ₁  -ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(n-1)*ExactReal(n-2)*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ end
                                    end
                                    if n-2 >= 0 && m-2 >= 0
                                        if i == 10 && j == 2  Z₁_V₄_ᵧ₁ = Z₁_V₄_ᵧ₁  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ end
                                        if i == 10 && j == 4  Z₁_V₄_ᵧ₃ = Z₁_V₄_ᵧ₃  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ end
                                    end
                                    if m-3 >=0
                                        if i == 10 && j == 4  Z₁_V₄_ᵧ₃ = Z₁_V₄_ᵧ₃  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ExactReal(m-1)*ExactReal(m-2)*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ end
                                    end
                                    if n-2 >= 0 && m-2>= 0 && ℓ >= 2
                                        if i == 10 && j == 1   Z₁_V₄_α = Z₁_V₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₂_int[n,m,ℓ] * (ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ + ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ ) end
                                    elseif n-2 >= 0 && ℓ >= 2
                                        if i == 10 && j == 1   Z₁_V₄_α = Z₁_V₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₂_int[n,m,ℓ] * ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1)*v₁ end
                                    elseif m-2 >= 0 && ℓ >= 2
                                        if i == 10 && j == 1   Z₁_V₄_α = Z₁_V₄_α - ExactReal(ℓ-1)*α_int^ExactReal(ℓ-2)*ω_int*M₂_int[n,m,ℓ] * ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2)*v₃ end
                                    end
                                end
                            end
                        end
                    end

                    if i == 5 && j == 2  Z₁_mat[i,j][k]  =  Λ_tails*norm(Z₁_Γ₄_ᵧ₁, 𝑿ᵧ ) end
                    if i == 5 && j == 4  Z₁_mat[i,j][k]  =  Λ_tails*norm(Z₁_Γ₄_ᵧ₃, 𝑿ᵧ ) end
                    if i == 5 && j == 24  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_Γ₄_ω, 𝑿ᵧ ) end
                    if i == 5 && j == 1  Z₁_mat[i,j][k]   =  Λ_tails*norm(Z₁_Γ₄_α, 𝑿ᵧ ) end

                    if i == 10 && j == 2  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₄_ᵧ₁, 𝑿ᵥ ) end
                    if i == 10 && j == 4  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₄_ᵧ₃, 𝑿ᵥ ) end
                    if i == 10 && j == 24  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₄_ω, 𝑿ᵥ )  end
                    if i == 10 && j == 1   Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₄_α, 𝑿ᵥ ) end
                    if i == 10 && j == 7  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_V₄_ᵥ₁, 𝑿ᵥ ) end
                    if i == 10 && j == 9  Z₁_mat[i,j][k]  = Λ_tails*norm(Z₁_V₄_ᵥ₃, 𝑿ᵥ ) end

                    if i == 14 && j == 11  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₄_w₁, 𝑿_w ) end
                    if i == 14 && j == 13  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₄_w₃, 𝑿_w ) end
                    if i == 14 && j == 24  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₄_ω, 𝑿_w )  end
                    if i == 14 && j == 1  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_W₄_α, 𝑿_w ) end

                    if i == 29 && j == 26  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₄_ₚ₁, 𝑿ₚ ) end
                    if i == 29 && j == 28  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₄_ₚ₃, 𝑿ₚ ) end
                    if i == 29 && j == 1  Z₁_mat[i,j][k] = Λ_tails*norm(Z₁_P₄_α, 𝑿ₚ ) end

                    if i == 22 && j == 20  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₄_ₐ₁, 𝑿ₐ )  end
                    if i == 22 && j == 22  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₄_ₐ₃, 𝑿ₐ ) end
                    if i == 22 && j == 1  Z₁_mat[i,j][k] = Λ_tails*abs(L)/ExactReal(4)*norm_T*norm(Z₁_G₄_α, 𝑿ₐ ) end
                    if i == 22 && j == 15  Z₁_mat[i,j][k] = Λ_tails/ExactReal(4)*norm_T*norm(Z₁_G₄_L, 𝑿ₐ ) end
                end
            end
            Z₁_[i,j] = Z₁_[i,j] + norm.(grid2cheb(Z₁_mat[i,j], DN[i,j]) , 1)[1]
        end
    end

    for i = 1:29
        #println("i = $(i)")
        if i in [1;6;15;16;17;18;19;24;25]
            local space_test_out = ℓ∞()
        elseif i in [ 2 ;3;4;5]
            local space_test_out =  ℓ¹(GeometricWeight(νᵧ))
        elseif i in [ 7;8;9;10]
            local space_test_out =  ℓ¹(GeometricWeight(νᵥ))
        elseif i in [ 11;12;13;14]
            local space_test_out =  ℓ¹(GeometricWeight(ν_w))
        elseif i in [ 20;21;22;23]
            local space_test_out =  ℓ¹(GeometricWeight(νₐ))
        elseif i in [ 26;27;28;29]
            local space_test_out =  ℓ¹(GeometricWeight(νₚ))
        end
        for k = 1:29
            #println("i = $(i) k = $(k)")
            if k in [1;6;15;16;17;18;23;24;25]
                local space_test_in = ℓ∞()
            elseif k in [ 2 ;3;4;5]
                local space_test_in =  ℓ¹(GeometricWeight(νᵧ))
            elseif k in [ 7;8;9;10]
                local space_test_in =  ℓ¹(GeometricWeight(νᵥ))
            elseif k in [ 11;12;13;14]
                local space_test_in =  ℓ¹(GeometricWeight(ν_w))
            elseif k in [ 19;20;21;22]
                local space_test_in =  ℓ¹(GeometricWeight(νₐ))
            elseif k in [ 26;27;28;29]
                local space_test_in =  ℓ¹(GeometricWeight(νₚ))
            end

            Aᴺ  =  LinearOperator(space_F_int[k] , space_X_int[i] , norm.(Aᵢⱼ_cheb[i,k] , 1))
            norm_Aᴺ = opnorm(Aᴺ, space_test_in, space_test_out)
            for j = 1:29
                #println("i = $(i) and j = $(j) and k = $(k)")
                if k ∈ [15;16;17;18] && j == k+5
                    Z₁_[i,j] = Z₁_[i,j] + norm_Aᴺ / ( νₐ^ExactReal(minimum(order(space_X_int_ext[j])+1)))
                end
                if k ∈ [15;16;17;18] && j == k+11
                    Z₁_[i,j] = Z₁_[i,j] + norm_Aᴺ * maximum( [ (norm_z/νₚ)^ExactReal.( order(space_X_int_ext[j])[1] + 1 ) , (norm_z/νₚ)^ExactReal.( order(space_X_int_ext[j])[2] + 1 ) ]) 
                end
                if k ∈ [19;20;21;22] && j == k+1
                    Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,:][:,1]) , space_test_out )/ ( νₐ^ExactReal(minimum(order(space_X_int_ext[j])+1) ))
                end
                if k ∈ [19;20;21;22] && j == k-8
                    Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,:][:,1]) , space_test_out ) * maximum( [ (norm_y/ν_w)^ExactReal.( order(space_X_int_ext[j])[1] + 1 ) , (ExactReal(0.95)/ν_w)^ExactReal.( order(space_X_int_ext[j])[2] +1 ) ])
                end
            end
        end
    end

    #=
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


        for j in [ 11;12;13;14;26;27;28;29;20;21;22;23] 
                if j == 11
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,19] , 1)
                elseif j == 12
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,20] , 1)
                elseif j == 13
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,21] , 1)
                elseif j == 14
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,22] , 1)
                elseif j == 26
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,15] , 1)
                elseif j == 27
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,16] , 1)
                elseif j == 28
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,17] , 1)
                elseif j == 29
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,18] , 1)
                elseif j == 20
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,19] , 1)
                elseif j == 21
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,20] , 1)
                elseif j == 22
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,21] , 1)
                elseif j == 23
                    Aᴺ  =  norm.(Aᵢⱼ_cheb[i,22] , 1)
                end

                if j == 11 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_y/ν_w)^( order_w₁[1] + 1 ) , (interval(0.95)/ν_w)^( order_w₁[2] +1 ) ]) end
                if j == 12 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_y/ν_w)^( order_w₂[1] + 1 ) , (interval(0.95)/ν_w)^( order_w₂[2] +1 ) ]) end
                if j == 13 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_y/ν_w)^( order_w₃[1] + 1 ) , (interval(0.95)/ν_w)^( order_w₃[2] +1 ) ]) end
                if j == 14 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_y/ν_w)^( order_w₄[1] + 1 ) , (interval(0.95)/ν_w)^( order_w₄[2] +1 ) ]) end
                
                if j == 26 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_z/νₚ)^( order_p₁[1] + 1 ) , (norm_z/νₚ)^( order_p₁[2] + 1 ) ]) end
                if j == 27 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_z/νₚ)^( order_p₂[1] + 1 ) , (norm_z/νₚ)^( order_p₂[2] + 1 ) ]) end
                if j == 28 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_z/νₚ)^( order_p₃[1] + 1 ) , (norm_z/νₚ)^( order_p₃[2] + 1 ) ]) end
                if j == 29 Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) .* maximum( [ (norm_z/νₚ)^( order_p₄[1] + 1 ) , (norm_z/νₚ)^( order_p₄[2] + 1 ) ]) end
                
                if j == 20 Z₁_[i,j] = Z₁_[i,j] + 2*norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) ./ ( νₐ^minimum(order_a₁)   ) end
                if j == 21 Z₁_[i,j] = Z₁_[i,j] + 2*norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) ./ ( νₐ^minimum(order_a₂)   ) end
                if j == 22 Z₁_[i,j] = Z₁_[i,j] + 2*norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) ./ ( νₐ^minimum(order_a₃)   ) end
                if j == 23 Z₁_[i,j] = Z₁_[i,j] + 2*norm( Sequence( space_X_int[i], Aᴺ[:,1]) , space_test_in ) ./ ( νₐ^minimum(order_a₄)   ) end
        end
    end

=#
    return Z₁_
end

function weight_x!(Z₁_ )
    vect_temp = abs.(real(eigvecs(sup.(Z₁_))[:,end]))
    weight_x = interval.(1 .* ( 1 ./   vect_temp ))
    return weight_x
end

function DF_ext_fft!(DF_ext_fft, X_cheb, η_cheb, DN_ftt, space_X_int, space_X_int_ext ,space_F_int, κ, star, lengthᵥ_int, M₁_int, M₂_int )
    iter = ProgressBar(1:29)
    set_description(iter, "    - ")
    for i in iter
        for j = 1:29
            if DN_ftt[i,j] > 1
                DF_ext_fft[i,j] =  [zeros(Complex{eltype(lengthᵥ_int)}, space_X_int_ext[j] , space_F_int[i] ) for _ ∈ 1:DN_ftt[i,j]]
                X_temp = cheb2grid(X_cheb, DN_ftt[i,j])
                η_temp = cheb2grid(η_cheb, DN_ftt[i,j])
                for ℓ = 1:DN_ftt[i,j]
                    DF_ext_fft[i,j][ℓ] = DFᵢⱼ(DF_ext_fft[i,j][ℓ],  Sequence(space_X_int , X_temp[ℓ]), κ, star, lengthᵥ_int, M₁_int, M₂_int ,  Sequence(space_X_int ,η_temp[ℓ])  , i, j )
                end
            end
        end
    end
    return DF_ext_fft
end


function DF_ext_fft_dft!(DF_ext_fft, X_cheb, η_cheb, DN, space_X_int, space_X_int_ext ,space_F_int, κ, star, lengthᵥ_int, M₁_int, M₂_int )
    for i = 1:29
        for j = 1:29
            if DN[i,j] > 1
                DF_ext_fft[i,j] = map(A -> interval.(A), [zeros(ComplexF64, space_X_int_ext[j] , space_F_int[i] ) for _ ∈ 1:(DN[i,j]+1)])
                X_temp = cheb2grid(X_cheb, DN[i,j])
                η_temp = cheb2grid(η_cheb, DN[i,j])
                for ℓ = 1:DN[i,j]+1
                    DF_ext_fft[i,j][ℓ] = DFᵢⱼ(DF_ext_fft[i,j][ℓ],  Sequence(space_X_int , X_temp[ℓ]), κ, star, lengthᵥ_int, M₁_int, M₂_int ,  Sequence(space_X_int ,η_temp[ℓ])  , i, j )
                end
            end
        end
    end
    return DF_ext_fft
end

function cheb2grid_LinearOperator(Aᵢⱼ_cheb, DN_ftt, i, j, k,space_F_int,space_X_int)
    test = cheb2grid(Aᵢⱼ_cheb[i,k], DN_ftt[k , j])
    return map(A -> LinearOperator( space_F_int[k], space_X_int[i], A), [ test[o]  for o ∈ 1:DN_ftt[k , j]])
end

function map_matrix(ADF_ext_fft_temp, DN_ftt, k, j)
    return map(A -> A[:,:], [ADF_ext_fft_temp[o] for o ∈ 1:DN_ftt[k , j]])
end
function map_matrix2(ADF_ext_fft_temp, DN_fttₖⱼ)
    return map(A -> A[:,:], [ADF_ext_fft_temp[o] for o ∈ 1:DN_fttₖⱼ])
end

function B!(B, Aᵢⱼ_cheb, DF_ext_fft,space_F_int, DN, DN_ftt,space_X_int)
    iter = ProgressBar(1:29)
    set_description(iter, "    ")
    for i in iter
        for j = 1:29
            ADF_ext_cheb = []
            for k = 1:29
                if DN_ftt[k,j] > 1
                    if isempty(ADF_ext_cheb )
                        ADF_ext_cheb = grid2cheb(map_matrix( cheb2grid_LinearOperator(Aᵢⱼ_cheb, DN_ftt, i, j, k,space_F_int,space_X_int).* DF_ext_fft[k , j] , DN_ftt, k, j) , DN[k,j] )
                        
                    else
                        ADF_ext_cheb = ADF_ext_cheb + grid2cheb(map_matrix( cheb2grid_LinearOperator(Aᵢⱼ_cheb, DN_ftt, i, j, k,space_F_int,space_X_int).* DF_ext_fft[k , j] , DN_ftt, k, j) , DN[k,j] )
                    end
                end
            end
            if i == j
                B[i,j] = abs.(UniformScaling(interval(1))  - LinearOperator( space_X_int[j], space_X_int[i], norm.(ADF_ext_cheb, 1) ) )
                #B[i,j] = abs.(UniformScaling(1)  - LinearOperator( space_X_int[j], space_X_int[i], norm.(ADF_ext_cheb, 1) ) )
            else
                B[i,j] =  LinearOperator( space_X_int[j], space_X_int[i], norm.(ADF_ext_cheb, 1) ) 
            end
        end
        # DF_ext_fft[:,j] .= nothing
        GC.gc()
    end
    return B
end

function B6!(B, Aᵢⱼ_cheb, DF_ext_fft,space_F_int, DN, DN_ftt,space_X_int, r, norm_A_mid, νᵧ , νᵥ, ν_w, νₐ, νₚ)
    checker = Matrix{Vector}(undef,29,29)
    for i = 1:29
        for k = 1:29
            checker[i,k] =  filter(x -> x ≠ 1, unique(DN_ftt[k,:]))
        end
    end
    iter = ProgressBar(1:29)
    set_description(iter, "    ")
    for i in iter
        for j = 1:29
            ADF_ext_cheb = []
            for k = 1:29
                if i in [1;6;15;16;17;18;19;24;25]
                    local space_test_in2 = ℓ∞()
                elseif i in [ 2 ;3;4;5]
                    local space_test_in2 =  ℓ¹(GeometricWeight(νᵧ))
                elseif i in [ 7;8;9;10]
                    local space_test_in2 =  ℓ¹(GeometricWeight(νᵥ))
                elseif i in [ 11;12;13;14]
                    local space_test_in2 =  ℓ¹(GeometricWeight(ν_w))
                elseif i in [ 20;21;22;23]
                    local space_test_in2 =  ℓ¹(GeometricWeight(νₐ))
                elseif i in [ 26;27;28;29]
                    local space_test_in2 =  ℓ¹(GeometricWeight(νₚ))
                end
                if k in [1;6;15;16;17;18;23;24;25]
                    local space_test_out2 = ℓ∞()
                elseif k in [ 2 ;3;4;5]
                    local space_test_out2 =  ℓ¹(GeometricWeight(νᵧ))
                elseif k in [ 7;8;9;10]
                    local space_test_out2 =  ℓ¹(GeometricWeight(νᵥ))
                elseif k in [ 11;12;13;14]
                    local space_test_out2 =  ℓ¹(GeometricWeight(ν_w))
                elseif k in [ 19;20;21;22]
                    local space_test_out2 =  ℓ¹(GeometricWeight(νₐ))
                elseif k in [ 26;27;28;29]
                    local space_test_out2 =  ℓ¹(GeometricWeight(νₚ))
                end
            
                if DN_ftt[k,j] > 1
                    A_temp = cheb2grid(Aᵢⱼ_cheb[i,k],  DN_ftt[k,j])
                    if DN_ftt[k,j] ∈ checker[i,k]
                        r[i,k] = max( r[i,k], interval.( maximum( maximum.( map.( X -> radius.(X), A_temp)))).* interval(-1,1))
                        A_temp = map(X -> mid.(X),A_temp)
                        norm_A_mid[i,k] = max( norm_A_mid[i,k], opnorm( LinearOperator( space_F_int[k],space_X_int[i] ,norm.(grid2cheb( map( X -> interval.(X), A_temp), DN_ftt[k,j]), 1)),space_test_out2,space_test_in2))
                        filter!(e->e≠DN_ftt[k,j],checker[i,k])
                    else
                        A_temp = map(X -> mid.(X),A_temp)
                    end
                    
                    if isempty(ADF_ext_cheb )
                        ADF_ext_cheb = grid2cheb( map( X -> interval(X), A_temp .* DF_ext_fft[k , j] ) , DN[k,j] )
                    else
                        ADF_ext_cheb = ADF_ext_cheb + grid2cheb( map( X -> interval(X), A_temp .* DF_ext_fft[k , j] ) , DN[k,j] )
                    end
                end
            end
            if i == j
                B[i,j] = LinearOperator(space_X_int[j] , space_X_int[i], norm.((UniformScaling(interval(1))  - LinearOperator( space_X_int[j], space_X_int[i], ADF_ext_cheb ))[:,:] , 1) ) 
            else
                B[i,j] =  LinearOperator( space_X_int[j] , space_X_int[i], norm.(LinearOperator( space_X_int[j], space_X_int[i], ADF_ext_cheb )[:,:] , 1)) 
            end
        end
        #Aᵢⱼ_cheb[i,:] .= nothing
        GC.gc()
    end
    return B 
end


function B7!(B, Aᵢⱼ_cheb, DF_ext_fft,space_F_int, DN, DN_ftt,space_X_int, r, norm_A_mid, νᵧ , νᵥ, ν_w, νₐ, νₚ)
    checker = Matrix{Vector}(undef,29,29)
    for i = 1:29
        for k = 1:29
            checker[i,k] =  filter(x -> x ≠ 1, unique(DN_ftt[k,:]))
        end
    end
    iter = ProgressBar(1:29)
    set_description(iter, "    ")
    for i in iter
        ADF_ext_cheb = Vector{Any}(undef,29)
        for k = 1:29
            if i in [1;6;15;16;17;18;19;24;25]
                local space_test_in2 = ℓ∞()
            elseif i in [ 2 ;3;4;5]
                local space_test_in2 =  ℓ¹(GeometricWeight(νᵧ))
            elseif i in [ 7;8;9;10]
                local space_test_in2 =  ℓ¹(GeometricWeight(νᵥ))
            elseif i in [ 11;12;13;14]
                local space_test_in2 =  ℓ¹(GeometricWeight(ν_w))
            elseif i in [ 20;21;22;23]
                local space_test_in2 =  ℓ¹(GeometricWeight(νₐ))
            elseif i in [ 26;27;28;29]
                local space_test_in2 =  ℓ¹(GeometricWeight(νₚ))
            end
            if k in [1;6;15;16;17;18;23;24;25]
                local space_test_out2 = ℓ∞()
            elseif k in [ 2 ;3;4;5]
                local space_test_out2 =  ℓ¹(GeometricWeight(νᵧ))
            elseif k in [ 7;8;9;10]
                local space_test_out2 =  ℓ¹(GeometricWeight(νᵥ))
            elseif k in [ 11;12;13;14]
                local space_test_out2 =  ℓ¹(GeometricWeight(ν_w))
            elseif k in [ 19;20;21;22]
                local space_test_out2 =  ℓ¹(GeometricWeight(νₐ))
            elseif k in [ 26;27;28;29]
                local space_test_out2 =  ℓ¹(GeometricWeight(νₚ))
            end
            N_check = length( checker[i,k] )
            for ℓ = 1:N_check
                A_temp = cheb2grid(Aᵢⱼ_cheb[i,k],  checker[i,k][ℓ])
                r[i,k] = max( r[i,k], interval.( maximum( maximum.( map.( X -> radius.(X), A_temp)))).* interval(-1,1))
                A_temp = map(X -> mid.(X),A_temp)
                norm_A_mid[i,k] = max( norm_A_mid[i,k], opnorm( LinearOperator( space_F_int[k],space_X_int[i] ,norm.(grid2cheb( map( X -> interval.(X), A_temp), checker[i,k][ℓ]), 1)),space_test_out2,space_test_in2))
                for j ∈ findall(x->x==checker[i,k][ℓ], DN_ftt[k,:])
                    if isassigned( ADF_ext_cheb , j)
                        ADF_ext_cheb[j] = ADF_ext_cheb[j] + grid2cheb( map( X -> interval(X), A_temp .* DF_ext_fft[k , j] ) , DN[k,j] )
                    else
                        ADF_ext_cheb[j] = grid2cheb( map( X -> interval(X), A_temp .* DF_ext_fft[k , j] ) , DN[k,j] )
                    end
                end
            end
        end
        for j = 1:29
            if i == j
                B[i,j] = LinearOperator(space_X_int[j] , space_X_int[i], norm.((UniformScaling(interval(1))  - LinearOperator( space_X_int[j], space_X_int[i], ADF_ext_cheb[j] ))[:,:] , 1) ) 
            else
                B[i,j] =  LinearOperator( space_X_int[j] , space_X_int[i], norm.(LinearOperator( space_X_int[j], space_X_int[i], ADF_ext_cheb[j] )[:,:] , 1)) 
            end
        end
        #Aᵢⱼ_cheb[i,:] .= nothing
        GC.gc()
    end
    return B 
end

function A_!(A_, Aᵢⱼ_cheb,space_F_int,space_X_int,νᵧ,νᵥ,ν_w,νₐ,νₚ)
    for i = 1:29
        for j = 1:29
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
            if j in [1;6;15;16;17;18;23;24;25]
                local space_test_out = ℓ∞()
            elseif j in [ 2 ;3;4;5]
                local space_test_out =  ℓ¹(GeometricWeight(νᵧ))
            elseif j in [ 7;8;9;10]
                local space_test_out =  ℓ¹(GeometricWeight(νᵥ))
            elseif j in [ 11;12;13;14]
                local space_test_out =  ℓ¹(GeometricWeight(ν_w))
            elseif j in [ 19;20;21;22]
                local space_test_out =  ℓ¹(GeometricWeight(νₐ))
            elseif j in [ 26;27;28;29]
                local space_test_out =  ℓ¹(GeometricWeight(νₚ))
            end
            Aᴺ = LinearOperator(  space_F_int[j] , space_X_int[i] , norm.(Aᵢⱼ_cheb[i,j],1) ) 
            A_[i, j ]  = A_[i, j ] + opnorm( Aᴺ , space_test_out , space_test_in )
        end
    end
    return A_
end

function finding_poly(a , b , r_star ,νₚ)

    init_ = a*r_star/νₚ
    initₙ = a*r_star*(a*r_star + 2*b)/(νₚ^2)
    count = 3
    while sup(init_) < inf(initₙ) && count <= 100
        init_ = copy(initₙ)
        initₙ = sup( ((Polynomial( [ b ; a ] , :r)^count - b^count)(r_star))/νₚ^count )

        count = count + 1
    end
    if count == 100
        error("finding_poly did not converge in less than 100 itterations ")
    end
    


end


function DN_ftt_index!(DN_ftt)
    DN_fft_index = []
    for i = 1:29
        for j = 1:29
            if DN_ftt[i,j] > 1 && DN_ftt[i,j] ∉ DN_fft_index
                DN_fft_index = [DN_fft_index ; DN_ftt[i,j] ] 
            end
        end
    end
    return DN_fft_index
end
function DN_index!(DN_ftt)
    DN_fft_index = []
    for i = 1:29
        for j = 1:29
            if DN_ftt[i,j] > 0 && DN_ftt[i,j] ∉ DN_fft_index
                DN_fft_index = [DN_fft_index ; DN_ftt[i,j] ] 
            end
        end
    end
    return DN_fft_index
end


function Ψ_cheb!(Ψ, a, NN ,M, ν)
    for k = 0:M
        Ψ_temp₁ = interval(0)
        Ψ_temp₂ = interval(0)
        if (k-NN) <= -(M+1)
            for j = (k-NN):-(M+1)
                Ψ_temp₁ = max( Ψ_temp₁, abs(a[abs.(k-j)])./ (ν ^ExactReal(abs(j))) )
            end
        end
        if (M+1) <= (k+NN)
            for j = (M+1):(k+NN)
                Ψ_temp₂ = max( Ψ_temp₂,  abs(a[abs.(k-j)])./ (ν ^ExactReal(abs(j))) )
            end
        end
        if eltype(ν) == Interval{Float64}
            Ψ[k] = interval.( max(  sup(Ψ_temp₁) ,  sup(Ψ_temp₂)  ))
        else
            Ψ[k] = max(  sup(Ψ_temp₁) ,  sup(Ψ_temp₂)  )
        end
    end
    return Ψ
end
function Ψ_Fourier!(Ψ, a, NN ,M, ν)
    for k = -M:M
        Ψ_temp₁ = interval(0)
        Ψ_temp₂ = interval(0)
        if (k-NN) <= -(M+1)
            for j = (k-NN):-(M+1)
                Ψ_temp₁ = max( Ψ_temp₁, abs(a[abs.(k-j)])./ (ν ^abs(j)) )
            end
        end
        if (M+1) <= (k+NN)
            for j = (M+1):(k+NN)
                Ψ_temp₂ = max( Ψ_temp₂,  abs(a[abs.(k-j)])./ (ν ^abs(j)) )
            end
        end
        Ψ[k] = max( Ψ_temp₁ , Ψ_temp₂  )
    end
    return Ψ
end

function Z₁_body!(Z₁_, X_cheb, Aᵢⱼ_cheb, DN, DN_ftt, M₁_int, M₂_int , νᵧ,  νᵥ, νₐ,  space_F ,space_X_int, ν_w, νₚ)
    norm_T = ExactReal(1)/νₐ + ExactReal(2)*νₐ
    zᵢⱼ_mat = Matrix{Vector}(undef, 29 , 29)
    Ψᵢⱼ_mat = Matrix{Vector}(undef, 29 , 29)
    iter = ProgressBar(1:29)
    set_description(iter, "    - ")
    for i in iter
        for j = 1:29
            zᵢⱼ_mat[i,j] = Vector{Sequence}(undef, DN_ftt[i,j] ) 
            Ψᵢⱼ_mat[i,j] = Vector{Any}(undef, DN_ftt[i,j] ) 
            if DN_ftt[i,j] > 1 && i ∈ [3;5;8;10;20;22]
                X_temp = cheb2grid( X_cheb   , DN_ftt[i,j]  )
                for k = 1:DN_ftt[i,j]
                     local X_ = Sequence( space_X_int , X_temp[k]  )
                     local α_int = real(X_[1])
                     local x = x_remove_complex!(component(X_,2:29))
                     local γ,λ,v,w,L,θ₀,σ,a,ω_int,eigenpairs,p = x2var!(x)
                     local γ₁,γ₂,γ₃,γ₄ = eachcomponent(γ)
                     local γ₁_ext,γ₂_ext,γ₃_ext,γ₄_ext = γ_extened_orientable!(γ)
                     local v₁,v₂,v₃,v₄ = eachcomponent(v)
                     local w₁, w₂, w₃, w₄ = eachcomponent(w)
                     local p₁, p₂, p₃, p₄ = eachcomponent(p)
                     local a₁, a₂, a₃, a₄ = eachcomponent(a)    
                    if i == 3 
                        if j == 2
                            zᵢⱼ_mat[i,j][k] = zero(γ₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),  space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0 
                                            if n >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵧ )[:]
                        elseif j == 4
                            zᵢⱼ_mat[i,j][k] = zero(γ₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0 
                                            if m >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k]  = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ]*ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵧ )[:]
                        end
                    elseif i == 5
                        if j == 2
                            zᵢⱼ_mat[i,j][k] = zero(γ₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0 
                                            if n >= 2 && (n-2 + m -1 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(n-1)*γ₁^(n-2)*γ₃^(m-1)
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵧ )[:]
                        elseif j == 4
                            zᵢⱼ_mat[i,j][k] = zero(γ₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0 
                                            if m >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k]  = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ]*ExactReal(m-1)*γ₁^(n-1)*γ₃^(m-2)
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵧ )[:]
                        end
                    elseif i == 8
                        if j == 7
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if n-2 >= 0 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ] *  ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) 
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        elseif j == 9
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if m-2 >= 0 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k]  = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₁_int[n,m,ℓ] *  ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) 
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        elseif j == 2
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if n-3 >= 0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  -ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ( ExactReal((n-1)*(n-2))*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ ) 
                                            end
                                            if n-2 >= 0 && m-2 >= 0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ ) 
                                            end
                
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]

                        elseif j == 4
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if n-2 >= 0 && m-2 >= 0 && (n + m -4) > 0
                                                Ψᵢⱼ_mat[i,j][k] = Ψᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ )
                                            end
                                            if m-3 >=0 && (n + m -4) > 0
                                                Ψᵢⱼ_mat[i,j][k] = Ψᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₁_int[n,m,ℓ] * ( ExactReal((m-1)*(m-2))*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ ) 
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        end
                    elseif i == 10
                        if j == 7
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if n-2 >= 0 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ] *  ExactReal(n-1)*γ₁_ext^(n-2)*γ₃_ext^(m-1) 
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        elseif j == 9
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if m-2 >= 0 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k]  = zᵢⱼ_mat[i,j][k]  - α_int^ExactReal(ℓ-1)*ω_int*M₂_int[n,m,ℓ] *  ExactReal(m-1)*γ₁_ext^(n-1)*γ₃_ext^(m-2) 
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        elseif j == 2
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if n-3 >= 0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  -ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ( ExactReal((n-1)*(n-2))*γ₁_ext^(n-3)*γ₃_ext^(m-1)*v₁ ) 
                                            end
                                            if n-2 >= 0 && m-2 >= 0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₃ ) 
                                            end
                
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        elseif j == 4
                            zᵢⱼ_mat[i,j][k] = zero(v₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if n-2 >= 0 && m-2 >= 0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ( ExactReal((m-1)*(n-1))*γ₁_ext^(n-2)*γ₃_ext^(m-2)*v₁ ) 
                                            end
                                            if m-3 >=0 && (n + m -4) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k]  - ω_int*α_int^ExactReal(ℓ-1)*M₂_int[n,m,ℓ] * ( ExactReal((m-1)*(m-2))*γ₁_ext^(n-1)*γ₃_ext^(m-3)*v₃ )
                                            end
                                        end
                                    end
                                end
                            end
                            Ψᵢⱼ_mat[i,j][k] = Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νᵥ )[:]
                        end
                    elseif i == 20
                        if j == 20
                            zᵢⱼ_mat[i,j][k] = zero(a₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if n >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] + α_int^ExactReal(ℓ-1)*L/ExactReal(4)*M₁_int[n,m,ℓ] * ExactReal(n-1)*a₁^(n-2)*a₃^(m-1)
                                            end
                                        end
                                    end
                                end
                            end  
                            Ψᵢⱼ_mat[i,j][k] = norm_T.*Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νₐ )[:]              
                        elseif j == 22
                            zᵢⱼ_mat[i,j][k] = zero(a₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₁_int,1)
                                for m in axes(M₁_int,2)
                                    for ℓ in axes(M₁_int,3)
                                        if M₁_int[n,m,ℓ] != 0
                                            if m >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] + α_int^ExactReal(ℓ-1)*L/ExactReal(4)*M₁_int[n,m,ℓ] * ExactReal(m-1)*a₁^(n-1)*a₃^(m-2)
                                            end
                                        end
                                    end
                                end
                            end  
                            Ψᵢⱼ_mat[i,j][k] = norm_T.*Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νₐ )[:]              
                        end
                    elseif i == 22
                        if j == 20
                            zᵢⱼ_mat[i,j][k] = zero(a₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros(eltype(α_int), space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if n >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] + α_int^ExactReal(ℓ-1)*L/ExactReal(4)*M₂_int[n,m,ℓ] * ExactReal(n-1)*a₁^(n-2)*a₃^(m-1)
                                            end
                                        end
                                    end
                                end
                            end  
                            Ψᵢⱼ_mat[i,j][k] = norm_T.*Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νₐ )[:]              
                        elseif j == 22
                            zᵢⱼ_mat[i,j][k] = zero(a₁) 
                            Ψᵢⱼ_mat[i,j][k] = zeros( eltype(α_int),space_F[i] )
                            for n in axes(M₂_int,1)
                                for m in axes(M₂_int,2)
                                    for ℓ in axes(M₂_int,3)
                                        if M₂_int[n,m,ℓ] != 0
                                            if m >= 2 && (n + m -3 ) > 0
                                                zᵢⱼ_mat[i,j][k] = zᵢⱼ_mat[i,j][k] + α_int^ExactReal(ℓ-1)*L/ExactReal(4)*M₂_int[n,m,ℓ] * ExactReal(m-1)*a₁^(n-1)*a₃^(m-2)
                                            end
                                        end
                                    end
                                end
                            end  
                            Ψᵢⱼ_mat[i,j][k] = norm_T.*Ψ_cheb!(Ψᵢⱼ_mat[i,j][k], zᵢⱼ_mat[i,j][k] , order(zᵢⱼ_mat[i,j][k]) ,order(space_F[i]), νₐ )[:]              
                        end
                    end
                end
            end
        end
    end


    AΨ_fft = Matrix{Any}(undef, 29,29)
    # q=5 and j=2  at 0.5
    for i = 1:29
        for j = 1:29
            for q = 1:29
                if DN_ftt[q,j] > 1 
                    if (q ∈ [3;5]  && j ∈ [2;4]) || (q ∈ [20;22]  && j ∈ [20;22]) || (q ∈ [8;10]  && j ∈ [2;4;7;9])
                        if isassigned(AΨ_fft,i,j)   
                            #AΨ_fft[i,j] = AΨ_fft[i,j] +  grid2cheb( dot.( abs_all(cheb2grid( Aᵢⱼ_cheb[i,q]   , DN_ftt[q,j]  ), DN_ftt[q,j]),Ψᵢⱼ_mat[q,j]) ,DN[q,j] )
                            AΨ_fft[i,j] = AΨ_fft[i,j] +  grid2cheb( dot2vec_of_mat( abs_all(cheb2grid( Aᵢⱼ_cheb[i,q] , DN_ftt[q,j]  ),DN_ftt[q,j]),  Ψᵢⱼ_mat[q,j])    , DN[q,j])
                        else
                            AΨ_fft[i,j] = grid2cheb( dot2vec_of_mat( abs_all(cheb2grid( Aᵢⱼ_cheb[i,q] , DN_ftt[q,j]  ),DN_ftt[q,j]), Ψᵢⱼ_mat[q,j])  , DN[q,j])
                        end
                    end
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
        for j = 1:29
            if isassigned(AΨ_fft,i,j) 
                Z₁_[i,j] = Z₁_[i,j] + norm( Sequence( space_X_int[i] , norm.( AΨ_fft[i,j] , 1 )[:] ) , space_test_in ) 
            end
        end
    end
    
    return Z₁_ 
end

function abs_all( A , NNN )
    A_abs = copy(A)
    for i = 1:NNN
        A_abs[i] = abs.(A[i])
    end
    return A_abs
end

function dot2vec_of_mat(A,B)
    if size(A[1],2) ==1
        C = dot.(A,B) 
        D = Vector{Matrix}(undef, length(C) )
        for i = 1:length(D)
            D[i] =  reshape( [C[i] ], 1,1)
        end
    else
        D = A .* transpose.(reshape.(B, 1,length(B[1])))
    end
    return D
 end

function test_A_cheb_error!(X_cheb, η_cheb , N, N_fft , space_X, space_F, κ, star, lengthᵥ, M₁, M₂ )
    X_temp = cheb2grid( X_cheb , N_fft )
    η_temp = cheb2grid( η_cheb , N_fft )
    A_fft = Vector{Matrix}(undef, N_fft )
    for i = 1:N_fft
        A_fft[i] = inv(DF_all!(zeros( ComplexF64 , space_X, space_F ) , Sequence( space_X , mid.( X_temp[i]) ), κ, star, lengthᵥ, M₁, M₂ , Sequence( space_X , mid.( η_temp[i]) ) ))[:,:]
    end
    A_cheb = grid2cheb( A_fft , N)
    A_fft_back = cheb2grid(A_cheb , N_fft)
    Δ = maximum(abs.( A_fft[1] -  A_fft_back[1] ) )
    return Δ
end

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

function cheb_remove_im!(X_cheb, X_grid)
    N_vec = Vector{Any}(undef, 29)
    for i = 1:29
        if i in [1;6;15;16;17;18;19;24]
            N_vec[i] = 1 
        elseif i in [25]
            N_vec[i] = 10 
        elseif i in [2;4]
            N_vec[i] = order(component(X_grid[1] ,i ))+1
        elseif i in [3;5]
            N_vec[i] = order(component(X_grid[1] ,i ))
        elseif i in [7;8;9;10]
            N_vec[i] = 2*order(component(X_grid[1] ,i ))+1
        elseif i in [11;12;13;14]
            N_vec[i] = (2*order(component(X_grid[1] ,i ))[1]+1)*(order(component(X_grid[1] ,i ))[2]+1)
        elseif i in [ 20;21;22;23]
            N_vec[i] = order(component(X_grid[1] ,i ))+1
        elseif i in [26;27;28;29]
            N_vec[i] = (order(component(X_grid[1] ,i ))[1]+1)*(order(component(X_grid[1] ,i ))[2]+1)
        end
    end
    M_vec = Vector{Any}(undef, 29)
    for i = 1:29
        M_vec[i] = sum(N_vec[1:i])
    end
    for k = 1:length(X_grid[1])
        if (  k <= M_vec[5] ) || ( k > M_vec[14]  && k <= M_vec[24]) 
            X_cheb[k] = real.(X_cheb[k])
        end
    end
    return X_cheb

end


function Vandermonde(t,y)
    V = zeros(eltype(t),length(t),length(t))
    V[:,1] .=  1
    for i = 2:length(t)
        V[:,i] = t.^(i-1)
    end
    return inv(V)*y
end

function Vandermonde_interval(t,y)
    V =  zeros(eltype(t),length(t),length(t))
    V[:,1] .=   1 
    for i = 2:length(t)
        V[:,i] = t.^(i-1) 
    end
    return interval.( inv(V) )*y
end


function grid2taylor(X_grid)
    N = length(X_grid)-1
    m =  length( X_grid[1])
    X_taylor = Vector{Any}(undef, m )
    t = range(-1,1,N+1) 
    for p = 1:m
        y = zeros(eltype( X_grid[1][p]) ,N+1,1)
        for ℓ  = 1:N+1
            y[ℓ] = X_grid[ℓ][p]
        end
        X_taylor[p] = Sequence( Taylor(N),Vandermonde(t,y)[:])
    end
    return X_taylor
end


function grid2taylor_int(X_grid)
    N = length(X_grid)-1
    m =  length( X_grid[1])
    X_taylor = Vector{Any}(undef, m )
    t =  range(-1,1,N+1) 
    for p = 1:m
        y = zeros(eltype( X_grid[1][p]) ,N+1,1)
        for ℓ  = 1:N+1
            y[ℓ] = X_grid[ℓ][p]
        end
        X_taylor[p] =  Sequence( Taylor(N),Vandermonde_interval(t,y)[:])
    end
    return X_taylor
end

function AF_taylor!( AF_grid)
    AFᵢ_taylor = Vector{Any}( undef, 29 )
    for i = 1:29
        AFᵢ_taylor[i] =  grid2taylor_int(AF_grid[i])
    end
    AF_taylor = copy(AFᵢ_taylor[1])
    for i = 2:29
        AF_taylor = AF_taylor + AFᵢ_taylor[i] 
    end
    return norm.(AF_taylor, 1)
end

function taylor2grid(X_Taylor, dN )
    N = length(X_Taylor[1])-1
    m =  length(X_Taylor)
    X_grid = Vector{Any}(undef, dN + 1  )
    t = range(-1,1, dN + 1 )
    for ℓ  = 1:dN+1
        X_grid[ℓ] = zeros(eltype(X_Taylor[1]),m,1)
        for p = 1:m
            X_grid[ℓ][p] = X_Taylor[p](t[ℓ])
        end
    end
    return X_grid
end


function taylor2grid_int(X_Taylor, dN )
    N = length(X_Taylor[1])-1
    m =  length(X_Taylor)
    X_grid = Vector{Any}(undef, dN + 1  )
    t = interval.( range(-1,1, dN + 1 ) )
    for ℓ  = 1:dN+1
        X_grid[ℓ] = zeros(eltype(X_Taylor[1]),m,1)
        for p = 1:m
            X_grid[ℓ][p] = X_Taylor[p](t[ℓ])
        end
    end
    return X_grid
end



function grid2taylor_int_A_col(A_grid)
    A_taylor = Vector{Any}(undef, 29 )
    for i = 1:29
        N = length(A_grid[i])-1
        t = range(-1,1,N+1) 
        m₁ =  size( A_grid[i][1],1)
        m₂ =  size( A_grid[i][1],2)
        A_taylor[i] = Matrix{Any}(undef, m₁, m₂ )
        for k = 1:m₁
            for j = 1:m₂
                y = zeros(eltype( A_grid[i][1][k,j]) ,N+1,1)
                for ℓ = 1:N+1
                    y[ℓ] = A_grid[i][ℓ][k,j]
                end
                A_taylor[i][k,j] = Sequence( Taylor(N),Vandermonde(t,y)[:])
            end
        end
        A_taylor[i] =  map(A -> interval.(A), A_taylor[i])
    end
    return A_taylor
end

function taylor2grid_int_A_col(A_Taylorᵢ, dNᵢ )
    t = interval.( range(-1,1, dNᵢ+1))
    A_grid_ext =  Vector{Any}(undef, dNᵢ+1)
    m₁ =  size( A_Taylorᵢ,1)
    m₂ =  size( A_Taylorᵢ,2)
    for ℓ = 1:dNᵢ+1
        A_grid_ext[ℓ] = zeros(eltype(A_Taylorᵢ[1]),m₁,m₂)
        for p₁ = 1:m₁
            for p₂ = 1:m₂
                A_grid_ext[ℓ][p₁,p₂] = A_Taylorᵢ[p₁,p₂](t[ℓ])
            end
        end
    end
    return A_grid_ext
end


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

function grid2cheb_dft(x_grid)
    Sequence(space(x_grid[1]), [idft!(getindex.(x_grid, i)) for i ∈ indices(space(x_grid[1]))])
end