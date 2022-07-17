using Random, LinearAlgebra

struct NortonSVIntegrator{TF,TG}
    dt::Float64
    η::Float64
    T::Float64
    γ::Float64
    F::TF # forcing profile
    G::TG # response profile
end
NortonSVIntegrator(dt::Float64, η::Float64, T::Float64, γ::Float64, F::Function, G::Function) = NortonSVIntegrator{typeof(F),typeof(G)}(dt, η, T, γ, F, G)

function Molly.simulate!(sys::System, sim::NortonSVIntegrator, n_steps; parallel::Bool=true, rng=Random.GLOBAL_RNG)
    force_hist=Float64[]

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt(1 - α^2)

    accels = accelerations(sys, neighbors; parallel=parallel)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot products
    FdotG = dot(F_y, G_y)
    GdotG = dot(G_y, G_y)

    #initialize state on constant response manifold
    λ = (sim.η-dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y
    
    run_loggers!(sys, neighbors, 0; parallel=parallel)
    for step_n = 1:n_steps
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_13 = (sim.η-dot(G_y,v_x))/FdotG#analytic expression for Lagrange multiplier
        v_x .+= λ_13 * F_y 

        #A step
        sys.coords .+= sys.velocities * sim.dt
        sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

        accels .= accelerations(sys, neighbors; parallel=parallel)

        F_y .= sim.F.(q_y)
        G_y .= sim.G.(q_y)

        FdotG = dot(F_y, G_y)
        GdotG = dot(G_y, G_y)

        (isnan(FdotG)) && return force_hist #abort if system NaNs out
        λ_12 = (sim.η-dot(G_y,v_x))/FdotG #correction term to reproject momenta on manifold
        v_x .+=λ_12 * F_y
        #λ_12/Δt is a good approximation to the term ∇R_q(q_t,p_t)⋅p_t/F(q_t)⋅G(q_t) forcing the dynamics to remain on the cotangent bundle
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_23 = (sim.η-dot(G_y,v_x))/FdotG
        v_x .+= λ_23 * F_y 
        
        #O_step
        velocities .= α*sys.velocities + σ * random_velocities(sys, sim.T; rng=rng)  #equilibrium fd solution
        λ_fd=(sim.η - dot(G_y,v_x))/FdotG #analytic expression for Lagrange multiplier
        v_x .+= λ_fd * F_y

        F_ham=(λ_13+λ_23)/sim.dt
        F_ou=sim.γ*sim.η/FdotG
        F_corr= λ_12/sim.dt 

        push!(force_hist, F_ham+F_corr+F_ou)

        neighbors = find_neighbors(sys, sys.neighbor_finder,neighbors,step_n; parallel=parallel)
        run_loggers!(sys, neighbors, step_n; parallel=parallel)
    end
    return force_hist
end