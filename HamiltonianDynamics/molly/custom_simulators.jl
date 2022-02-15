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


function Molly.simulate!(sys::System{D,S,false},
    sim::SymplecticEulerA,
    n_steps::Integer; parallel::Bool = true) where {D,S}

    if any(inter -> !inter.nl_only, values(sys.general_inters))
        neighbors_all = Molly.all_neighbors(length(sys))
    else
        neighbors_all = nothing
    end
    neighbors = find_neighbors(sys, sys.neighbor_finder, parallel = parallel)

    accels = zero(sys.velocities)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        for i = 1:length(sys)
            sys.coords[i] += sys.velocities[i] * sim.dt
            sys.coords[i] = wrap_coords.(sys.coords[i], sys.box_size)
        end
        accels = Molly.remove_molar.(accelerations(sys, neighbors, parallel = parallel))

        for i = 1:length(sys)
            sys.velocities[i] += accels[i] * sim.dt
        end

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel = parallel)
        end
    end
    return sys
end


function Molly.simulate!(sys::System{D,S,false},
    sim::SymplecticEulerB,
    n_steps::Integer;
    parallel::Bool = true) where {D,S}

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

        for i = 1:length(sys)
            sys.velocities[i] += accels[i] * sim.dt
            sys.coords[i] += sys.velocities[i] * sim.dt
            sys.coords[i] = wrap_coords.(sys.coords[i], sys.box_size)
        end

        neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, parallel = parallel)
    end
    return sys
end

struct ExplicitEuler{T,C}
    dt::T
    coupling::C
end

ExplicitEuler(; dt, coupling = NoCoupling()) = ExplicitEuler(dt, coupling)

function Molly.simulate!(sys::System{D,S,false},
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
        sys.coords += sys.velocities .* sim.dt
        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
        sys.velocities += accels * sim.dt

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n)
        end
    end
    return sys
end
