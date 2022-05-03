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
    dW = zero(sys.coords)
    candidate_coords = zero(sys.coords)
    accels = accelerations(sys, neighbors; parallel=parallel)
    accels_tilde = zero(accels)
    V = potential_energy(sys, neighbors)
    V_tilde = zero(V)
    σ = sqrt(2 * sim.dt)

    reverse_deviation = zero(sys.coords)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)
        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        ## propose candidate position according to Euler-Maruyama scheme ##
        @. candidate_coords = sys.coords + sim.β * accels * sim.dt + σ * dW
        @. reverse_deviation = sys.coords - candidate_coords #aperiodic part of the reverse displacement

        ## temporarily swap candidate position in the system's field to compute candidate potential and gradient ##
        candidate_coords.=wrap_coords_vec.(candidate_coords,(sys.box_size,))#enforce boundary conditions
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
        exp_α=exp(α)
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

mutable struct ImprovedMALA
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
end

function ImprovedMALA(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    ImprovedMALA(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function Molly.simulate!(sys::System{D}, sim::ImprovedMALA, n_steps::Integer; parallel::Bool=true) where {D} end