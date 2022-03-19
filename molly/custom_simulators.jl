struct SymplecticEulerA{T,C}
    dt::T
    coupling::C
end

SymplecticEulerA(; dt, coupling=NoCoupling()) = SymplecticEulerA(dt, coupling)

struct SymplecticEulerB{T,C}
    dt::T
    coupling::C
end

SymplecticEulerB(; dt, coupling=NoCoupling()) = SymplecticEulerB(dt, coupling)


function Molly.simulate!(sys::System{D},
    sim::SymplecticEulerA,
    n_steps::Integer; parallel::Bool=true) where {D}

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel=parallel)

    accels = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        @. sys.coords += sys.velocities * sim.dt

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel=parallel))

        @. sys.velocities += accels * sim.dt

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel=parallel)
        end
    end
    return sys
end


function Molly.simulate!(sys::System{D},
    sim::SymplecticEulerB,
    n_steps::Integer;
    parallel::Bool=true) where {D}

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel=parallel)
    accels = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel=parallel))

        @. sys.velocities += accels * sim.dt
        @. sys.coords += sys.velocities * sim.dt

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel=parallel)
    end
    return sys
end

struct ExplicitEuler{T,C}
    dt::T
    coupling::C
end

ExplicitEuler(; dt, coupling=NoCoupling()) = ExplicitEuler(dt, coupling)

function Molly.simulate!(sys::System{D,false},
    sim::ExplicitEuler,
    n_steps::Integer;
    parallel::Bool=true) where {D,S}

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)
        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel=parallel))
        @. sys.coords += sys.velocities * sim.dt
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
        @. sys.velocities += accels * sim.dt

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n)
        end
    end
    return sys
end

struct LangevinBAO
    dt::Real
    γ::Real
    β::Real
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBAO(; dt, γ, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))

    β = inv(T) #todo work with units, i.e. kb !=1
    LangevinBAO(dt, γ, β, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBAO,
    n_steps::Integer;
    parallel::Bool=true) where {D}

    M_inv = inv.(mass.(sys.atoms))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)
    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        @. sys.velocities += accels_t * sim.dt#B

        @. sys.coords += sys.velocities * sim.dt#A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.velocities = α * sys.velocities + σ * dW#O

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)

        accels_t = accelerations(sys, neighbors; parallel=parallel)
    end

end

struct LangevinBAOA
    dt::Real
    γ::Real
    β::Real
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBAOA(; dt, γ, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))

    β = inv(T) #todo work with units, i.e. kb !=1
    LangevinBAOA(dt, γ, β, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBAOA,
    n_steps::Integer;
    parallel::Bool=true) where {D}

    M_inv = inv.(mass.(sys.atoms))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)
        @. sys.velocities += accels_t * sim.dt#B
        @. sys.coords += sys.velocities * sim.dt / 2#A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.velocities = α * sys.velocities + σ * dW#O
        @. sys.coords += sys.velocities * sim.dt / 2#A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)

        accels_t = accelerations(sys, neighbors; parallel=parallel)
    end

end


struct LangevinBABO
    dt::Real
    γ::Real
    β::Real
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBABO(; dt, γ, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))

    β = inv(T) #todo work with units, i.e. kb !=1
    LangevinBABO(dt, γ, β, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBABO,
    n_steps::Integer;
    parallel::Bool=true) where {D}

    M_inv = inv.(mass.(sys.atoms))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)
    @showprogress for step_n in 1:n_steps

      run_loggers!(sys, neighbors, step_n)

        @. sys.velocities += accels_t * sim.dt / 2 #B 


        @. sys.coords += sys.velocities * sim.dt #A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        accels_t .= accelerations(sys, neighbors; parallel=parallel)


        @. sys.velocities += accels_t * sim.dt / 2 #B 


        dW .= SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.velocities = α * sys.velocities + σ * dW #O

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
    end

end


struct LangevinBAOAB
    dt::Real
    γ::Real
    β::Real
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBAOAB(; dt, γ, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))

    β = ustrip(inv(T)) #todo work with units, i.e. kb !=1
    LangevinBAOAB(dt, γ, β, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBAOAB,
    n_steps::Integer;
    parallel::Bool=true) where {D}

    M_inv = inv.(ustrip.(mass.((sys.atoms))))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)


        @. sys.velocities += accels_t * sim.dt / 2 #B

        @. sys.coords += sys.velocities * sim.dt / 2 #A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.velocities = α * sys.velocities + σ * dW #O

        @. sys.coords += sys.velocities * sim.dt / 2 #A
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        accels_t = accelerations(sys, neighbors; parallel=parallel)

        @. sys.velocities += accels_t * sim.dt / 2 #B

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
    end

end

struct LangevinSplitting
    dt::Real
    γ::Real
    β::Real

    rseed::UInt32
    rng::AbstractRNG

    splitting::AbstractString
end

function LangevinSplitting(; dt, γ, T, splitting, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T)) #todo work with units, i.e. kb !=1
    @assert all(x ∈ "ABO" for x in splitting) "Invalid splitting descriptor: use only letters A, B and O."
    LangevinSplitting(dt, γ, β, rseed, rng, splitting)
end


function Molly.simulate!(sys::System{D}, sim::LangevinSplitting, n_steps::Integer, parallel::Bool=true) where {D}

    M_inv = inv.(ustrip.(mass.((sys.atoms))))

    α_eff = zero(M_inv)
    σ_eff = zero(M_inv)

    @. α_eff = exp(-sim.γ * sim.dt * M_inv / count('O', sim.splitting))
    @. σ_eff = sqrt(M_inv * (1 - α_eff^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.pairwise_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    effective_dts = [sim.dt / count(c, sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    #determine the need to recompute accelerations before B steps

    #first pass to determine if first B step needs to compute accelerations
    for op in sim.splitting
        if op == 'A'
            forces_known = false
        elseif op == 'B'
            forces_known = true
        end
    end

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    steps = []
    arguments = []

    for (j, op) in enumerate(sim.splitting)
        if op == 'A'
            push!(steps, A_step!)
            push!(arguments, (sys, effective_dts[j]))
        elseif op == 'B'
            push!(steps, B_step!)
            push!(arguments, (sys, effective_dts[j], accels_t, neighbors, force_computation_steps[j], parallel))
        elseif op == 'O'
            push!(steps, O_step!)
            push!(arguments, (sys, α_eff, σ_eff, sim.rng, dW))
        end
    end

    step_arg_pairs = zip(steps, arguments)

    for step_n = 1:n_steps

        run_loggers!(sys, neighbors, step_n)

        for (step!, args) = step_arg_pairs

            step!(args...)

        end

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
    end
end

function O_step!(s::System{D}, α_eff::Vector{T}, σ_eff::Vector{T}, rng::AbstractRNG, noise_vec::Vector{SVector{D,T}}) where {D,T}
    noise_vec .= SVector{D}.(eachrow(randn(rng, Float64, (length(sys), D))))
    @. s.velocities = α_eff * s.velocities + σ_eff * noise_vec
end

function A_step!(s::System, dt_eff::Real)
    @. s.coords += s.velocities * dt_eff
    s.coords .= wrap_coords_vec.(s.coords, (s.box_size,))
end

function B_step!(s::System{D}, dt_eff::Real, acceleration_vector::Vector{SVector{D,T}}, neighbors, compute_forces::Bool, parallel::Bool) where {D,T}
    compute_forces && (acceleration_vector .= accelerations(s, neighbors, parallel=parallel))

    @. s.velocities += dt_eff * acceleration_vector
end

mutable struct LangevinGHMC
    dt::Real
    γ::Real
    β::Real

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
end

function LangevinGHMC(; dt, γ, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))

    β = ustrip(inv(T)) #todo work with units, i.e. kb !=1
    LangevinGHMC(dt, γ, β, rseed, rng, 0, 0)
end




function Molly.simulate!(sys::System{D}, sim::LangevinGHMC, n_steps::Integer, parallel::Bool=true) where {D}
    M_inv = inv.(ustrip.(mass.((sys.atoms))))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    candidate_neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors; parallel=parallel)

    accels_t = zero(sys.velocities)
    accels_t_dt = zero(sys.velocities)

    dW = zero(sys.velocities)

    candidate_coords = zero(sys.coords)
    candidate_velocities = zero(sys.velocities)

    accels_t = accelerations(sys, neighbors; parallel=parallel)

    for i = 1:n_steps
        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
        @. sys.velocities = α * sys.velocities + σ * dW
        H = total_energy(sys, neighbors)
        @. candidate_coords = sys.coords + sys.velocities * sim.dt + accels_t * sim.dt^2 / 2
        candidate_coords = wrap_coords_vec.(candidate_coords, (sys.box_size,))

        sys.coords, candidate_coords = candidate_coords, sys.coords

        @. candidate_velocities = sys.velocities + (accels_t + accels_t_dt) * sim.dt / 2

        sys.velocities, candidate_velocities = candidate_velocities, sys.velocities
        candidate_neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
        H_tilde = total_energy(sys, candidate_neighbors)

        U = rand(sim.rng)

        if min(1, exp(sim.β * (H - H_tilde))) > U
            accels_t = accels_t_dt
            neighbors = candidate_neighbors
            sim.n_accepted += 1
        else
            sys.coords, candidate_coords = candidate_coords, sys.coords
            sys.velocities, candidate_velocities = candidate_velocities, sys.velocities
            @. sys.velocities = -sys.velocities
        end

        sim.n_total += 1
    end

end