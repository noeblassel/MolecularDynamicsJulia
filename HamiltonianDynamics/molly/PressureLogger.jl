struct PressureLoggerJL{T}
    n_steps::Int
    pressures::Vector{T}
    r_c::Real
    ρ::Real
end

PressureLoggerJL(T, n_steps::Integer, ρ) = PressureLoggerJL(n_steps, T[], 2.5, ρ)
PressureLoggerJL(n_steps::Integer, ρ) = PressureLoggerJL(n_steps, Float64[], 2.5, ρ)
PressureLoggerJL(n_steps::Integer, r_c::Real, ρ) = PressureLoggerJL(n_steps, Float64[], r_c, ρ)


function Molly.log_property!(logger::PressureLoggerJL, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        N = length(s)
        r_c2 = logger.r_c^2
        K = Molly.kinetic_energy_noconvert(s)

        W = 0
        for i = 1:N
            for j = 1:i-1
                r_ij2 = Molly.square_distance(i, j, s.coords, s.box_size)
                if r_ij2 < r_c2
                    W += Molly.force_divr_nocutoff(s.general_inters[1], r_ij2, inv(r_ij2), (1.0, 1.0)) * r_ij2 #force_divr_nocutoff(LJ,r2,1/r2,(σ,ϵ))=-LJ'(r)/r
                end
            end
        end

        P = logger.ρ * (K / 1.5N + W / 3N)
        push!(logger.pressures, P)
    end
end