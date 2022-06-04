export
    SymplecticEulerA,
    SymplecticEulerB,
    ExplicitEuler,
    LangevinBAO,
    LangevinBAOA,
    LangevinBABO,
    LangevinBAOAB,
    LangevinSplitting,
    LangevinGHMC,
    MALA,
    MALA_HMC,
    NortonTest,
    NortonHomogeneousSplitting,
    NortonOneDriftSplitting,
    NortonTwoDriftSplitting,
    NortonColorDriftSplitting,
    NortonShearViscosityTest

log_barker(α::Float64) = (α > 0) ? (-α - log(1 + exp(-α))) : (-log(1 + exp(α)))
sq_norm(v::Vector{SVector{3,Float64}}) = dot(v, v)

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

function LangevinSplitting(; dt, γ, T, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
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

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

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
    accels_tilde = zero(accels)

    H = total_energy(sys, neighbors)
    H_tilde = zero(H)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
        @. sys.velocities = α * sys.velocities + σ * dW

        @. candidate_velocities = sys.velocities + accels * sim.dt / 2

        @. candidate_coords = sys.coords + candidate_velocities * sim.dt
        candidate_coords = wrap_coords_vec.(candidate_coords, (sys.box_size,))
        sys.coords, candidate_coords = candidate_coords, sys.coords

        accels_tilde = accelerations(sys, neighbors; parallel=parallel)

        @. candidate_velocities += accels_tilde * sim.dt / 2
        sys.velocities, candidate_velocities = candidate_velocities, sys.velocities

        H_tilde = total_energy(sys, neighbors)

        U = rand(sim.rng)

        sim.n_total += 1

        if log(U) < -sim.β * ustrip(H_tilde - H)
            sim.n_accepted += 1
            H = H_tilde
            accels, accels_tilde = accels_tilde, accels
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, i; parallel=parallel)
        else
            sys.coords, candidate_coords = candidate_coords, sys.coords
            sys.velocities, candidate_velocities = candidate_velocities, sys.velocities
            @. sys.velocities = -sys.velocities
        end

    end

end

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

struct NortonTest #simple BAO test in the one drift diagonal mass case
    dt::Real
    γ::Real
    β::Real
    v::Real

    rseed::UInt32
    rng::AbstractRNG
end

function NortonTest(; dt, γ, T, v, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    β = inv(T)
    return NortonTest(dt, γ, β, v, rseed, rng)
end

function Molly.simulate!(sys::System{D}, sim::NortonTest, n_steps::Integer; parallel::Bool=true) where {D}
    sys.velocities[1] = SVector(sim.v, 0.0, 0.0) + sys.velocities[1] .* SVector(0.0, 1.0, 1.0) #initialize state
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt((1 - α^2) / sim.β)
    N = length(sys)
    for step_n = 1:n_steps
        run_loggers!(sys, neighbors, step_n)
        accels = accelerations(sys, neighbors; parallel=parallel)
        for i = 2:N
            sys.velocities[i] += sim.dt * accels[i]
        end
        sys.velocities[1] += sim.dt * (accels[1] .* SVector(0.0, 1.0, 1.0))
        sys.coords += sim.dt * sys.velocities
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        for i = 2:N
            sys.velocities[i] = α * sys.velocities[i] + σ * dW[i]
        end
        sys.velocities[1] = (sys.velocities[1] .* SVector(1.0, 0.0, 0.0)) + ((α * sys.velocities[1] + σ * dW[1]) .* SVector(0.0, 1.0, 1.0))
        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
    end
end

struct NortonHomogeneousSplitting
    dt::Real
    γ::Real
    β::Real
    v::Real

    F::Vector{SVector{3,Float64}}
    rseed::UInt32
    rng::AbstractRNG

    splitting::AbstractString
end

function NortonHomogeneousSplitting(; dt, γ, T, v, F, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    β = ustrip.(inv.(T))
    @assert (all(x ∈ "ABO" for x ∈ splitting) && all(x ∈ splitting for x ∈ "ABO")) "Invalid splitting descriptor: use only and all letters A, B and O."
    NortonHomogeneousSplitting(dt, γ, β, v, F, rseed, rng, splitting)
end


function OneDriftSplitting(; N, dt, γ, T, v, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    one_drift = OneDriftNEMD(N, 1.0).force_field
    return NortonHomogeneousSplitting(; N, dt, γ, T, v, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
end
function NortonColorDriftSplitting(; N, dt, γ, T, v, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    color_drift = ColorDriftNEMD(N, 1.0).force_field
    return NortonHomogeneousSplitting(dt=dt, γ=γ, T=T, v=v, F=color_drift, splitting=splitting, rseed=rseed, rng=rng)
end

function NortonTwoDriftSplitting(; N, dt, γ, T, v, splitting, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    two_drift = TwoDriftNEMD(N, 1.0).force_field
    return NortonHomogeneousSplitting(dt=dt, γ=γ, T=T, v=v, F=two_drift, splitting=splitting, rseed=rseed, rng=rng)
end


"""Splitting scheme for mobility Norton dynamics, in the case of identity mass matrix."""
function Molly.simulate!(sys::System{D}, sim::NortonHomogeneousSplitting, n_steps::Integer, parallel::Bool=true) where {D}

    O_count = count('O', sim.splitting)
    B_count = count('B', sim.splitting)
    A_count = count('A', sim.splitting)

    dt_O = sim.dt / O_count
    dt_A = sim.dt / B_count
    dt_B = sim.dt / A_count

    α_eff = exp(-sim.γ * dt_O)
    σ_eff = sqrt((1 - α_eff^2) / sim.β)

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    accels = accelerations(sys, neighbors; parallel=parallel)
    inv_sq_norm_F = inv(dot(sim.F, sim.F))
    P(v::Vector{SVector{D,Float64}}) = sim.F * (dot(sim.F, v)) * inv_sq_norm_F#projector onto F
    P_perp(v::Vector{SVector{D,Float64}}) = v - P(v)#projector orthogonal to F
    F_coord(v::Vector{SVector{D,Float64}}) = dot(sim.F, v) * inv_sq_norm_F
    sys.velocities = sim.v * sim.F + P_perp(sys.velocities) #project initial state onto constant response hyperplane

    function A_step_F!()
        sys.coords += dt_A * sys.velocities
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
    end

    function O_step_F!()
        G = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
        sys.velocities = sim.v * sim.F + P_perp(α_eff * sys.velocities + σ_eff * G)
    end

    function B_step_F!(compute_forces::Bool, parallel::Bool=true)
        compute_forces && (accels .= accelerations(sys, neighbors, parallel=parallel))
        sys.velocities += dt_B * P_perp(accels)
    end

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

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
            push!(steps, A_step_F!)
            push!(arguments, ())
        elseif op == 'B'
            push!(steps, B_step_F!)
            push!(arguments, (force_computation_steps[j], parallel))
        elseif op == 'O'
            push!(steps, O_step_F!)
            push!(arguments, ())
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

struct NortonShearViscosityTest
    dt::Real
    γ::Real
    β::Real
    v::Real#response intensity

    F::Function #forcing profile
    rseed::UInt32
    rng::AbstractRNG
end

function NortonShearViscosityTest(; dt::Real, γ::Real, T::Real, v::Real, F::Function, rseed=round(UInt32, time()), rng=MersenneTwister(rseed))
    β = inv(T)
    return NortonShearViscosityTest(dt, γ, β, v, F, rseed, rng)
end

function Molly.simulate!(sys::System{D}, sim::NortonShearViscosityTest, n_steps::Integer; parallel::Bool=true) where {D}
    N = length(sys)

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt((1 - α^2) / sim.β)

    for step_n = 1:n_steps
        run_loggers!(sys, neighbors, step_n)
        accels = accelerations(sys, neighbors; parallel=parallel)

        for i = 1:N #B step
            F_y = sim.F(sys.coords[i][2])
            sys.velocities[i] = SVector{D,Float64}(vcat(sim.v * F_y, sys.velocities[i][2:end] + (sim.dt / 2) * accels[i][2:end]))
        end

        sys.coords += sim.dt * sys.velocities #A step
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        accels = accelerations(sys, neighbors; parallel=parallel)

        for i = 1:N #B step
            F_y = sim.F(sys.coords[i][2])
            sys.velocities[i] = SVector{D,Float64}(vcat(sim.v * F_y, sys.velocities[i][2:end] + (sim.dt / 2) * accels[i][2:end]))
        end

        for i = 1:N#O step
            G = randn(sim.rng, Float64, D - 1)
            sys.velocities[i] = SVector{D,Float64}(vcat(first(sys.velocities[i]), α * sys.velocities[i][2:end] + σ * G))
        end
        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
        (isnan(dot(sys.velocities,sys.velocities))) && (println("EXPLOSION at step $step_n !!!");exit(1))
    end
end