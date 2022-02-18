using NBodySimulator, Base.Threads

@fastmath function pair_virial_lj(sr::SimulationResult, time::Real, parallel::Bool = true)

    lj_params = sr.simulation.system.potentials[:lennard_jones]
    σ6 = lj_params.σ^2
    ϵ = lj_params.ϵ
    pbc = sr.simulation.boundary_conditions

    N = length(res.simulation.system.bodies)
    q = reshape(get_position(sr, time), (N, 3))

    if parallel && nthreads() > 1

        Ws_thread = [0.0 for i in 1:nthreads()]

        @threads for i = 1:N
            for j = 1:i-1
                ri = q[i]
                rj = q[j]
                (_, _, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

                if r2 < lj_params.R2
                    rm6 = inv(r2^3)

                end

            end
        end

        return sum(Ws_thread)


    else
        W = 0
        for i = 1:N
            for j = 1:i-1
                ri = q[i]
                rj = q[j]
                (_, _, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

                if r2 < lj_params.R2

                end

            end
        end

    end

end