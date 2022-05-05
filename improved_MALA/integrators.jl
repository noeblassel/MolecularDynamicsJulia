using Molly


log_barker(α::Float64) = (α > 0) ? (-α - log(1 + exp(-α))) : (-log(1 + exp(α)))
sq_norm(v::Vector{SVector{3,Float64}}) = dot(v, v)

mutable struct MALA
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
    baker_abs_sum::Float64
end

function MALA(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    MALA(dt, β, is_metropolis, rseed, rng, 0, 0, 0.0)
end

function Molly.simulate!(sys::System{D}, sim::MALA, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    F = forces(sys, neighbors; parallel=parallel)#compute the gradient, assumes all atoms have mass 1
    V = potential_energy(sys, neighbors)
    σ = sqrt(2 * sim.dt)
    λ = (σ * sim.β) / 2
    coords_tmp = zero(sys.coords)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)

        coords_tmp .= sys.coords #record original coordinates

        G = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        ## propose candidate position according to Euler-Maruyama scheme ##
        sys.coords += sim.β * sim.dt * F + σ * G
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))#enforce boundary conditions

        V_tilde = potential_energy(sys, neighbors)#candidate potential
        F_tilde = forces(sys, neighbors; parallel=parallel)#candidate gradient

        ## compute minus log of metropolis ratio ##
        α = sim.β * (V_tilde - V) + sq_norm(G + λ * (F + F_tilde)) / 2 - sq_norm(G) / 2 # log metropolis ratio

        ## accept / reject step ##
        U = rand(sim.rng)
        accepted = (sim.is_metropolis) ? (log(U) < -α) : (log(U) < log_barker(α)) #Metropolis-Hastings  / Barker rule
        sim.n_total += 1
        (!sim.is_metropolis) && (sim.baker_abs_sum += abs(2 * exp(-α) / (1 + exp(-α)) - 1))

        if accepted
            sim.n_accepted += 1
            F .= F_tilde #reuse candidate gradient as new gradient
            V = V_tilde #reuse candidate potential as new potential

            neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
        else
            sys.coords .= coords_tmp #revert to original coordinates
        end
    end

end

mutable struct MALA_HMC
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
    baker_abs_sum::Float64
end

function MALA_HMC(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    MALA_HMC(dt, β, is_metropolis, rseed, rng, 0, 0, 0.0)
end

function Molly.simulate!(sys::System{D}, sim::MALA_HMC, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    V = potential_energy(sys, neighbors)

    h = sqrt(2 * sim.β * sim.dt) #effective timestep of the symplectic scheme
    σ = inv(sqrt(sim.β)) #noise scale

    sys.velocities = σ * SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D)))) #sample random momentum
    coords_tmp = zero(sys.coords)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)


        # keep track of original position and Hamiltonian
        coords_tmp .= sys.coords
        H = V + Molly.kinetic_energy_noconvert(sys)

        # position Verlet evolution
        sys.coords += h * sys.velocities / 2 #A
        sys.velocities += h * accelerations(sys, neighbors) #B (only place gradient needs to be computed)
        sys.coords += h * sys.velocities / 2 #A

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))


        #compute new potential
        V_tilde = potential_energy(sys, neighbors)
        α = sim.β * (V_tilde + Molly.kinetic_energy_noconvert(sys) - H) # log of metropolis ratio rewrites as difference of Hamiltonians

        ## accept / reject step ##
        U = rand(sim.rng)

        accepted = (sim.is_metropolis) ? (log(U) < -α) : (log(U) < log_barker(α))
        sim.n_total += 1
        (!sim.is_metropolis) && (sim.baker_abs_sum += abs(2 * exp(-α) / (1 + exp(-α)) - 1))

        if accepted
            sim.n_accepted += 1
            neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel) #recompute neighbors
            V = V_tilde #reuse candidate potential as new potential
        else
            sys.coords .= coords_tmp
        end

        #sample new momentum for next iteration
        sys.velocities = σ * SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
    end

end