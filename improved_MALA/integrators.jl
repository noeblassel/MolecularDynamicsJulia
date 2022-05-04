using Molly

mutable struct MALA
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
end

function MALA(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    MALA(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function Molly.simulate!(sys::System{D}, sim::MALA, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    accels = accelerations(sys, neighbors; parallel=parallel)
    V = potential_energy(sys, neighbors)
    σ = sqrt(2 * sim.dt)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)
        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        ## propose candidate position according to Euler-Maruyama scheme ##
        candidate_coords = sys.coords + sim.β * accels * sim.dt + σ * dW
        displacement = candidate_coords - sys.coords #keep track of the displacement
        candidate_coords = wrap_coords_vec.(candidate_coords, (sys.box_size,))#enforce boundary conditions

        ## temporarily swap candidate position in the system's field to compute candidate potential and gradient ##
        tmp_coords=copy(sys.coords)
        sys.coords.=candidate_coords

        V_tilde = potential_energy(sys, neighbors)
        accels_tilde = accelerations(sys, neighbors; parallel=parallel)

        ## compute log of metropolis ratio ##
        α = sim.β * (V_tilde - V) # potential part
        +dot(displacement + sim.β * sim.dt * accels_tilde, displacement + sim.β * sim.dt * accels_tilde) / (4 * sim.dt) #reverse transition term
        -dot(displacement - sim.β * sim.dt * accels, displacement - sim.β * sim.dt * accels) / (4 * sim.dt) # forward transition term

        ## accept / reject step ##
        U = rand(sim.rng)
        accepted = (sim.is_metropolis) ? (log(U) < α) : (exp_α = exp(-α);(U < (exp_α / (1 + exp_α)))) #Metropolis-Hastings  / Barker rule
        sim.n_total += 1

        if accepted
            sim.n_accepted += 1
            neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel) #recompute neighbors
            accels, accels_tilde = accels_tilde, accels #reuse candidate gradient as new gradient
            V = V_tilde #reuse candidate potential as new potential
        else
            sys.coords.=tmp_coords #swap back original coordinates
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
end

function MALA_HMC(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    MALA_HMC(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function Molly.simulate!(sys::System{D}, sim::MALA_HMC, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    dW = zero(sys.coords)
    candidate_coords = zero(sys.coords)
    accels = accelerations(sys, neighbors; parallel=parallel)
    accels_mid = zero(accels)
    accels_tilde = zero(accels)
    V = potential_energy(sys, neighbors)
    V_tilde = zero(V)
    h = sqrt(2 * sim.β * sim.dt)
    σ = inv(sqrt(sim.β))
    displacement = zero(sys.coords)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)
        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        #put proposal coordinates in the system to use Molly's potential_energy and acceleration functions
        candidate_coords, sys.coords = sys.coords, candidate_coords

        ## propose candidate position according to position Verlet scheme from random momentum ##
        sys.velocities = σ * dW #random momentum
        sys.coords += h * sys.velocities / 2 #A
        accels_mid = accelerations(sys, neighbors)
        sys.velocities += h * accels_mid #B
        sys.coords += h * sys.velocities / 2 #A
        accels_tilde = accelerations(sys, neighbors)

        @. displacement = candidate_coords - sys.coords #aperiodic part of the reverse displacement

        ## temporarily swap candidate position in the system's field to compute candidate potential and gradient ##
        candidate_coords .= wrap_coords_vec.(candidate_coords, (sys.box_size,))#enforce boundary conditions
        sys.coords, candidate_coords = candidate_coords, sys.coords

        V_tilde = potential_energy(sys, neighbors)
        accels_tilde .= accelerations(sys, neighbors; parallel=parallel)

        ## update reverse displacement with gradient part ##
        @. reverse_deviation -= sim.β * sim.dt * accels_tilde

        ## compute log of metropolis ratio ##
        α = -sim.β * (V_tilde - V) # potential part
        +dot(dW, dW) / 2 # forward transition term
        -dot(reverse_deviation, reverse_deviation) / (4 * sim.dt) #reverse transition term

        ## accept / reject step ##
        U = rand(sim.rng)
        exp_α = exp(α)
        accepted = (sim.is_metropolis) ? (log(U) < α) : (U < (exp_α / (1 + exp_α)))
        sim.n_total += 1

        if accepted
            sim.n_accepted += 1
            neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel) #recompute neighbors
            accels, accels_tilde = accels_tilde, accels #reuse candidate gradient as new gradient
            V = V_tilde #reuse candidate potential as new potential
        else
            sys.coords, candidate_coords = candidate_coords, sys.coords #swap back original coordinates
        end
    end

end