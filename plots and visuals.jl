# ================================
# Functions for plotting
# ================================

#Plots
function plot_hetero_orbit(γ,a,p,w)
    Nₜ = 1000
    t = range(0,2*pi,Nₜ)
    tₐ = range(-1,1,Nₜ)
    γ_graph = zeros(Nₜ,4)
    v_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    for i = 1:Nₜ
        γ_graph[i,:] = γ(t[i])[:]
        a_graph[i,:] = a(tₐ[i])[:]
    end
    N_w = 200
    t_w = range( -1,1,N_w)
    r = range( 0,1,N_w) 
    theta_p = range( 0,2*pi,N_w)
    theta_w = range( 0,4*pi,N_w)
    w_graph = zeros(N_w,N_w,4)
    p_graph = zeros(N_w,N_w,4)

    for i = 1:N_w
        for j = 1:N_w
            w_graph[i,j,:] = real( w(theta_w[i],t_w[j]) )[:]
            σ₁_grpah = r[j]*cos( theta_p[i] ) + 1im*r[j]*sin( theta_p[i] )
            σ₂_grpah = r[j]*cos( theta_p[i] ) - 1im*r[j]*sin( theta_p[i] )
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
        end
    end
     P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.75,2), ylimits = (-1.5,1.5), zlimits = (-1.5,1.25)) #SH Snaking

    # P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.5,1.5), ylimits = (-1,1), zlimits = (-1,1))#SH Isolas
    
    # P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.2,0.4), ylimits = (-0.1,0.1), zlimits = (-1,0.6))# GS Snaking

    plot3d!( w_graph[:,:,1][:],w_graph[:,:,2][:],w_graph[:,:,3][:],linecolor = "red",alpha = 0.025 , title = "")
    annotate!((0, 0.25, 2, text("α = $(1)", :bottom, 18)))
    plot3d!( p_graph[:,:,1][:],p_graph[:,:,2][:],p_graph[:,:,3][:],linecolor = "cyan",alpha = 0.25 )
    plot3d!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:], linewidth=2, color =  "mediumorchid4")
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")

    return P
end
function plot_hetero_orbit2(γ,a,p,w)
    Nₜ = 1000
    t = range(0,2*pi,Nₜ)
    tₐ = range(-1,1,Nₜ)
    γ_graph = zeros(Nₜ,4)
    v_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    for i = 1:Nₜ
        γ_graph[i,:] = γ(t[i])[:]
        a_graph[i,:] = a(tₐ[i])[:]
    end
    N_w = 200
    t_w = range( -1,1,N_w)
    r = range( 0,1,N_w) 
    theta_p = range( 0,2*pi,N_w)
    theta_w = range( 0,4*pi,N_w)
    w_graph = zeros(N_w,N_w,4)
    p_graph = zeros(N_w,N_w,4)

    for i = 1:N_w
        for j = 1:N_w
            w_graph[i,j,:] = real( w(theta_w[i],t_w[j]) )[:]
            σ₁_grpah = r[j]*cos( theta_p[i] ) + 1im*r[j]*sin( theta_p[i] )
            σ₂_grpah = r[j]*cos( theta_p[i] ) - 1im*r[j]*sin( theta_p[i] )
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
        end
    end
    P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750) ,xlimits = (-0.2,0.7), ylimits = (-0.6,0.6), zlimits = (-0.1,0.1))
    plot3d!( w_graph[:,:,1][:],w_graph[:,:,2][:],w_graph[:,:,3][:],linecolor = "red",alpha = 0.25 )
    plot3d!( p_graph[:,:,1][:],p_graph[:,:,2][:],p_graph[:,:,3][:],linecolor = "cyan",alpha = 0.5 )
    plot3d!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:],linewidth=2,linecolor = "mediumorchid",alpha = 1 )
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")

    return P
end

function plot_hetero_orbit3(γ,a,p,w,α)
    Nₜ = 100
    t = range(0,2*pi,Nₜ)
    tₐ = range(-1,1,Nₜ)
    γ_graph = zeros(Nₜ,4)
    v_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    for i = 1:Nₜ
        γ_graph[i,:] = γ(t[i])[:]
        a_graph[i,:] = a(tₐ[i])[:]
    end
    N_w = 100
    t_w = range( -1,1,N_w)
    r = range( 0,1,N_w) 
    theta_p = range( 0,2*pi,N_w)
    theta_w = range( 0,4*pi,N_w)
    w_graph = zeros(N_w,N_w,4)
    p_graph = zeros(N_w,N_w,4)

    for i = 1:N_w
        for j = 1:N_w
            w_graph[i,j,:] = real( w(theta_w[i],t_w[j]) )[:]
            σ₁_grpah = r[j]*cos( theta_p[i] ) + 1im*r[j]*sin( theta_p[i] )
            σ₂_grpah = r[j]*cos( theta_p[i] ) - 1im*r[j]*sin( theta_p[i] )
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
        end
    end
    
    P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.75,2), ylimits = (-1.5,1.5), zlimits = (-1.5,1.25)) #SH Snaking
   
    t = range(0, 2π, length=200)
    s = range(-1, 1, length=100)
    r = range(0, 1, length=100)

    T = [ti for ti in t, si in s]
    S = [si for ti in t, si in s]
    R = [ri for ti in t, ri in r]

    x =  real( component(w,1).(T,S) )
    y = real( component(w,2).(T,S)  )
    z = real( component(w,3).(T,S)  )

    u =  real( component(p,1).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )
    v = real( component(p,2).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )
    uu = real( component(p,3).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )

   
    plot!( x, y, z, alpha = 0.05, lw = 2, color = :red , title = "" , fillcolor = :red)
    annotate!((0, 0.25, 2, text("α = $(round(α, digits = 4))", :bottom, 18)))
    plot!( u,v,uu, lw = 2,linecolor = "blue",alpha = 0.05 , fillrange=[u u] )
    plot!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:], linewidth=2, color =  "mediumorchid4")
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")

    return P
end

function plot_hetero_orbit4(γ,a,p,w,α)
    Nₜ = 1000
    t = range(0,2*pi,Nₜ)
    tₐ = range(-1,1,Nₜ)
    γ_graph = zeros(Nₜ,4)
    v_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    for i = 1:Nₜ
        γ_graph[i,:] = γ(t[i])[:]
        a_graph[i,:] = a(tₐ[i])[:]
    end
    N_w = 200
    t_w = range( -1,1,N_w)
    r = range( 0,1,N_w) 
    theta_p = range( 0,2*pi,N_w)
    theta_w = range( 0,4*pi,N_w)
    w_graph = zeros(N_w,N_w,4)
    p_graph = zeros(N_w,N_w,4)

    for i = 1:N_w
        for j = 1:N_w
            w_graph[i,j,:] = real( w(theta_w[i],t_w[j]) )[:]
            σ₁_grpah = r[j]*cos( theta_p[i] ) + 1im*r[j]*sin( theta_p[i] )
            σ₂_grpah = r[j]*cos( theta_p[i] ) - 1im*r[j]*sin( theta_p[i] )
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
        end
    end
    # P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.75,2), ylimits = (-1.5,1.5), zlimits = (-1.5,1.25)) #SH Snaking

     #P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.5,1.5), ylimits = (-1,1), zlimits = (-1,1))#SH Isolas
    
    P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.2,0.4), ylimits = (-0.1,0.1), zlimits = (-1,0.6))# GS Snaking

    t = range(0, 2π, length=200)
    s = range(-1, 1, length=100)
    r = range(0, 1, length=100)

    T = [ti for ti in t, si in s]
    S = [si for ti in t, si in s]
    R = [ri for ti in t, ri in r]

    x =  real( component(w,1).(T,S) )
    y = real( component(w,2).(T,S)  )
    z = real( component(w,3).(T,S)  )

    u =  real( component(p,1).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )
    v = real( component(p,2).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )
    uu = real( component(p,3).(R .* cos.(T) .+ 1im.* R .* sin.(T) , R .* cos.(T) .- 1im.* R .* sin.(T) )  )

   
    plot!( x, y, z, alpha = 0.025, lw = 1, color = :red , title = "" , fillcolor = :red)
    annotate!((0, 0.25, 2, text("α = $(round(α, digits = 4))", :bottom, 18)))
    plot!( u,v,uu, lw = 2,linecolor = "blue",alpha = 0.05 , fillrange=[u u] )
    plot!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:], linewidth=2, color =  "mediumorchid4")
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")

    return P
end

function plot_hetero_orbit5(γ,a,p,w,α)
    Nₜ = 1000
    t = range(0,2*pi,Nₜ)
    tₐ = range(-1,1,Nₜ)
    γ_graph = zeros(Nₜ,4)
    v_graph = zeros(Nₜ,4)
    a_graph = zeros(Nₜ,4)
    for i = 1:Nₜ
        γ_graph[i,:] = γ(t[i])[:]
        a_graph[i,:] = a(tₐ[i])[:]
    end
    N_w = 200
    t_w = range( -1,1,N_w)
    r = range( 0,1,N_w) 
    theta_p = range( 0,2*pi,N_w)
    theta_w = range( 0,4*pi,N_w)
    w_graph = zeros(N_w,N_w,4)
    p_graph = zeros(N_w,N_w,4)

    for i = 1:N_w
        for j = 1:N_w
            w_graph[i,j,:] = real( w(theta_w[i],t_w[j]) )[:]
            σ₁_grpah = r[j]*cos( theta_p[i] ) + 1im*r[j]*sin( theta_p[i] )
            σ₂_grpah = r[j]*cos( theta_p[i] ) - 1im*r[j]*sin( theta_p[i] )
            p_graph[i,j,:] = real( p(σ₁_grpah,σ₂_grpah) )[:]
        end
    end
    # P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.75,2), ylimits = (-1.5,1.5), zlimits = (-1.5,1.25)) #SH Snaking

     P  =  plot3d( [γ_graph[:,1]],[γ_graph[:,2]],[γ_graph[:,3]],linecolor = "red",linewidth=2, legend = false,xlabel = "u₁",ylabel = "u₂",zlabel = "u₃",camera = (60, 40),size=(750,750),xlimits = (-0.5,1.5), ylimits = (-1,1), zlimits = (-1,1))#SH Isolas
   
    t = range(0, π, length=200)
    s = range(-1, 1, length=100)
    tr = range(0, 2π, length=200)
    r = range(0, 1, length=100)

    T = [ti for ti in t, si in s]
    S = [si for ti in t, si in s]
    Tr = [tri for tri in tr, ri in r]
    R = [ri for tri in tr, ri in r]

    x =  real( component(w,1).(T,S) )
    y = real( component(w,2).(T,S)  )
    z = real( component(w,3).(T,S)  )

    u =  real( component(p,1).(R .* cos.(Tr) .+ 1im.* R .* sin.(Tr) , R .* cos.(Tr) .- 1im.* R .* sin.(Tr) )  )
    v = real( component(p,2).(R .* cos.(Tr) .+ 1im.* R .* sin.(Tr) , R .* cos.(Tr) .- 1im.* R .* sin.(Tr) )  )
    uu = real( component(p,3).(R .* cos.(Tr) .+ 1im.* R .* sin.(Tr) , R .* cos.(Tr) .- 1im.* R .* sin.(Tr) )  )

 
    plot!( x, y, z, alpha = 0.05, lw = 2, color = :red , title = "" , fillcolor = :red)
    annotate!((0, 0.25, 2, text("α = $(round(α, digits = 4))", :bottom, 18)))
    plot!( u,v,uu, lw = 1,linecolor = "blue",alpha = 0.05 , fillrange=[u u] )
    plot!( a_graph[:,1][:],a_graph[:,2][:],a_graph[:,3][:], linewidth=2, color =  "mediumorchid4")
    scatter!([0], [0],[0], color = "blue", label = "", markersize = 5 ,markerstrokecolor ="blue")

    return P
end