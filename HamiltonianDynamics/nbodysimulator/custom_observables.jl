LennardJonesParameters
@fastmath function pair_virial_lj(sr::NBodySimulator.SimulationResult, time::Real, parallel::Bool = true)

    lj_params = sr.simulation.system.potentials[:lennard_jones]
    σm6 = inv(lj_params.σ2^3)
    ϵ = lj_params.ϵ
    R2=lj_params.R2
    pbc = sr.simulation.boundary_conditions

    N = length(sr.simulation.system.bodies)
    q = get_position(sr, time)

    if parallel && nthreads() > 1

        W_threads = [0.0 for i in 1:nthreads()]

        @threads for i = 1:N
            for j = 1:i-1
                ri = q[:,i]
                rj = q[:,j]
                (_, _, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

                if r2 < R2
                    rm6 = inv(r2^3)
                    wij=24ϵ*σm6*(2σm6*rm6^2-rm6)
                    W_threads[threadid()]+=wij
                end

            end
        end

        return sum(W_threads)
    else
        W = 0
        for i = 1:N
            for j = 1:i-1
                ri = @view q[:,i]
                rj = @view q[:,j]
                (_, _, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

                if r2<R2
                    rm6 = inv(r2^3)
                    wij=24ϵ*σm6*(2σm6*rm6^2-rm6)
                    W+=wij
                end

            end
        end
        return W
    end

end