
function sim_lennard_jones_fluid(N_per_side, ρ, T, Δt, steps; species = nothing, reduced_units = true)
    if reduced_units
        ρʳ = ρ
        Tʳ = T
    else
        @assert species !== nothing "cannot compute reduced units for unspecified species"
        ρʳ = get_reduced_density(species, ρ)
        Tʳ = get_reduced_temperature(species, T)
    end

    N = N_per_dim^3
    L = (N / ρʳ)^(1 // 3)
    bodies = generate_bodies_in_cell_nodes(N, 1.0, sqrt(Tʳ), L)
    potential_parameters = LennardJonesParameters(1.0, 1.0, 2.5)

    system = PotentialNBodySystem(bodies, Dict(:lennard_jones => potential_parameters))
    boundary_conditions = CubicPeriodicBoundaryConditions(L)
    sim = NBodySimulation(system, (0.0, steps * Δt), boundary_conditions, NBodySimulator.NullThermostat(), 1.0)
    @time result = run_simulation(sim, VelocityVerlet(), dt = Δt)

    return result
end