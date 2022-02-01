function simulate!(sys::System{D, false},
                    sim::SymplecticEuler,
                    n_steps::Integer;
                    parallel::Bool=true) where {D, S}
    # See https://www.saylor.org/site/wp-content/uploads/2011/06/MA221-6.1.pdf for
    #   integration algorithm - used shorter second version
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    accels_t = accelerations(sys, neighbors; parallel=parallel)
    accels_t_dt = zero(accels_t)

    @showprogress for step_n in 1:n_steps
        run_loggers!(sys, neighbors, step_n)

        # Update coordinates
        for i in 1:length(sys)
            sys.coords[i] += sys.velocities[i] * sim.dt + remove_molar(accels_t[i]) * (sim.dt ^ 2) / 2
            sys.coords[i] = wrap_coords.(sys.coords[i], sys.box_size)
        end

        accels_t_dt = accelerations(sys, neighbors; parallel=parallel)

        # Update velocities
        for i in 1:length(sys)
            sys.velocities[i] += remove_molar(accels_t[i] + accels_t_dt[i]) * sim.dt / 2
        end

        apply_coupling!(sys, sim, sim.coupling)

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n; parallel=parallel)
            accels_t = accels_t_dt
        end
    end
    return sys
end