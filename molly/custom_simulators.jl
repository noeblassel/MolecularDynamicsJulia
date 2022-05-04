export 
    SymplecticEulerA,
    SymplecticEulerB,
    ExplicitEuler,
    LangevinBAO,
    LangevinBAOA,
    LangevinBABO,
    LangevinBAOAB,
    LangevinSplitting,
    LangevinGHMC

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

    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel=parallel)

    accels = zero(sys.velocities)

    for step_n in 1:n_steps
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

    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel=parallel)
    accels = zero(sys.velocities)

    for step_n in 1:n_steps
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

    neighbors = find_neighbors(sys, sys.neighbor_finder)

    for step_n in 1:n_steps
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

    M_inv = inv.(ustrip.(mass.(sys.atoms)))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)
    for step_n in 1:n_steps
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

    M_inv = inv.(ustrip.(mass.(sys.atoms)))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    for step_n in 1:n_steps
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

    M_inv = inv.(ustrip.(mass.(sys.atoms)))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)
    for step_n in 1:n_steps

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

    M_inv = inv.(ustrip.(mass.(sys.atoms)))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    for step_n in 1:n_steps
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

struct LangevinSplitting{T}
    dt::Real
    γ::Real
    β::T

    rseed::UInt32
    rng::AbstractRNG

    splitting::AbstractString
end

function LangevinSplitting(; dt, γ, T, splitting, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip.(inv.(T))
    @assert (all(x ∈ "ABO" for x ∈ splitting) && all(x ∈ splitting for x ∈ "ABO")) "Invalid splitting descriptor: use only and all letters A, B and O."
    LangevinSplitting{typeof(β)}(dt, γ, β, rseed, rng, splitting)
end


function Molly.simulate!(sys::System{D}, sim::LangevinSplitting, n_steps::Integer, parallel::Bool=true) where {D}

    M_inv = inv.(ustrip.(mass.(sys.atoms)))

    α_eff = zero(M_inv)
    σ_eff = zero(M_inv)

    @. α_eff = exp(-sim.γ * sim.dt * M_inv / count('O', sim.splitting))
    @. σ_eff = sqrt(M_inv * (1 - α_eff^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    accels_t = accelerations(sys, neighbors; parallel=parallel)
    dW = zero(sys.velocities)

    effective_dts = [sim.dt / count(c, sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$",sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

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

    d_start=Dates.now()

    for step_n = 1:n_steps

        run_loggers!(sys, neighbors, step_n)

        for (step!, args) = step_arg_pairs

            step!(args...)

        end

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
    end
end

function O_step!(s::System{D}, α_eff::V, σ_eff::V, rng::AbstractRNG, noise_vec::N) where {D,T,N<:AbstractVector{SVector{D,T}},V<:AbstractVector{T}}
    noise_vec .= SVector{D}.(eachrow(randn(rng, Float64, (length(s), D))))
    s.velocities = α_eff .* s.velocities + σ_eff .* noise_vec
end

function A_step!(s::System, dt_eff::Real)
    s.coords += s.velocities * dt_eff
    s.coords = wrap_coords_vec.(s.coords, (s.box_size,))
end

function B_step!(s::System{D}, dt_eff::Real, acceleration_vector::A, neighbors, compute_forces::Bool, parallel::Bool) where {D,T,A<:AbstractVector{SVector{D,T}}}
    compute_forces && (acceleration_vector .= accelerations(s, neighbors, parallel=parallel)) 
    s.velocities += dt_eff * acceleration_vector
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

function Molly.simulate!(sys::System{D}, sim::LangevinGHMC, n_steps::Integer; parallel::Bool=true) where {D}
    M_inv = inv.(ustrip.(mass.((sys.atoms))))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    dW = zero(sys.velocities)

    candidate_coords = zero(sys.coords)
    candidate_velocities = zero(sys.velocities)

    accels = accelerations(sys, neighbors; parallel=parallel)
    accels_tilde=zero(accels)
    
    H = total_energy(sys, neighbors)
    H_tilde=zero(H)
    
    for i=1:n_steps
        run_loggers!(sys, neighbors, i)

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
        @. sys.velocities = α * sys.velocities + σ * dW 

        @. candidate_velocities=sys.velocities + accels * sim.dt / 2

        @. candidate_coords = sys.coords + candidate_velocities * sim.dt
        candidate_coords = wrap_coords_vec.(candidate_coords, (sys.box_size,))
        sys.coords, candidate_coords = candidate_coords, sys.coords
        
        accels_tilde = accelerations(sys,neighbors; parallel=parallel)
        
        @. candidate_velocities += accels_tilde * sim.dt / 2
        sys.velocities, candidate_velocities = candidate_velocities, sys.velocities

        H_tilde = total_energy(sys, neighbors)

        U = rand(sim.rng)
        
        sim.n_total += 1

        if log(U) < -sim.β * ustrip(H_tilde-H)
            sim.n_accepted += 1
            H=H_tilde
            accels,accels_tilde=accels_tilde,accels
            neighbors=find_neighbors(sys, sys.neighbor_finder,neighbors,i; parallel=parallel)
        else
            sys.coords, candidate_coords = candidate_coords, sys.coords
            sys.velocities, candidate_velocities = candidate_velocities, sys.velocities
            @. sys.velocities = -sys.velocities 
        end

    end

end

