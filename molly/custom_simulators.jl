struct SymplecticEulerA{T,C}
    dt::T
    coupling::C
end

SymplecticEulerA(; dt, coupling = NoCoupling()) = SymplecticEulerA(dt, coupling)

struct SymplecticEulerB{T,C}
    dt::T
    coupling::C
end

SymplecticEulerB(; dt, coupling = NoCoupling()) = SymplecticEulerB(dt, coupling)


function Molly.simulate!(sys::System{D},
    sim::SymplecticEulerA,
    n_steps::Integer; parallel::Bool = true) where {D}

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel = parallel)

    accels = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        @. sys.coords += sys.velocities * sim.dt

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel = parallel))

        @. sys.velocities += accels * sim.dt

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel = parallel)
        end
    end
    return sys
end


function Molly.simulate!(sys::System{D},
    sim::SymplecticEulerB,
    n_steps::Integer;
    parallel::Bool = true) where {D}

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel = parallel)
    accels = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel = parallel))

        @. sys.velocities += accels * sim.dt
        @. sys.coords += sys.velocities * sim.dt

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel = parallel)
    end
    return sys
end

struct ExplicitEuler{T,C}
    dt::T
    coupling::C
end

ExplicitEuler(; dt, coupling = NoCoupling()) = ExplicitEuler(dt, coupling)

function Molly.simulate!(sys::System{D,false},
    sim::ExplicitEuler,
    n_steps::Integer;
    parallel::Bool = true) where {D,S}

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)
        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel = parallel))
        @. sys.coords += sys.velocities * sim.dt
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
        @. sys.velocities += accels * sim.dt

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n)
        end
    end
    return sys
end

struct LangevinBABO{C}
    dt::Real
    γ::Real
    β::Real
    coupling::C
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBABO(; dt, γ, T, coupling = NoCoupling(), rseed = UInt32(round(time())), rng = MersenneTwister(rseed))

    β = inv(T) #todo work with units, i.e. kb !=1
    LangevinBABO{typeof(coupling)}(dt, γ, β, coupling, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBABO,
    n_steps::Integer;
    parallel::Bool = true) where {D}

    M_inv = inv.(mass.(sys.atoms))

    α = zero(M_inv)
    σ = zero(M_inv)

    @. α = exp(-sim.γ * sim.dt * M_inv)
    @. σ = sqrt(M_inv * (1 - α^2) / sim.β) #noise on velocities, not momenta

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel = parallel)

    accels_t = accelerations(sys, neighbors; parallel = parallel)
    accels_t_dt = zero(accels_t)
    dW = zero(sys.velocities)
    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        @. sys.coords += sys.velocities * sim.dt + accels_t * sim.dt^2 / 2
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        accels_t_dt = accelerations(sys, neighbors; parallel = parallel)

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.velocities = α * (sys.velocities + (accels_t + accels_t_dt) * sim.dt / 2) + σ * dW

        accels_t = accels_t_dt
        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel = parallel)
    end

end


struct LangevinBAOAB{C}
    dt::Real
    γ::Real
    β::Real
    coupling::C
    rseed::UInt32
    rng::AbstractRNG
end


function LangevinBAOAB(; dt, γ, T, coupling = NoCoupling(), rseed = UInt32(round(time())), rng = MersenneTwister(rseed))

    β = ustrip(inv(T)) #todo work with units, i.e. kb !=1
    LangevinBAOAB{typeof(coupling)}(dt, γ, β, coupling, rseed, rng)
end

function Molly.simulate!(sys::System{D},
    sim::LangevinBAOAB,
    n_steps::Integer;
    parallel::Bool = true) where {D}

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

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel = parallel)

    accels_t = accelerations(sys, neighbors; parallel = parallel)
    accels_t_dt = zero(accels_t)
    dW = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))

        @. sys.coords += (1 + α) * (sys.velocities + accels_t * sim.dt / 2) * sim.dt / 2 + σ * dW * sim.dt / 2
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

        accels_t_dt = accelerations(sys, neighbors; parallel = parallel)

        @. sys.velocities = α * sys.velocities + (α * accels_t + accels_t_dt) * sim.dt / 2 + σ * dW

        accels_t = accels_t_dt
        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel = parallel)
    end

end

struct MALA{C}
    dt::Real
    γ::Real
    β::Real
    coupling::C
    rseed::UInt32
    rng::AbstractRNG
    α::Real
    σ::Real
end