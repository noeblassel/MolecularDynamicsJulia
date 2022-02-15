
function sim_lennard_jones_fluid(N_per_dim, ρ, T, Δt, steps,integrator)
    N = N_per_dim^3
    L = (N / ρ)^(1 // 3)
    bodies = generate_bodies_in_cell_nodes(N, 1.0, sqrt(T), L)
    potential_parameters = LennardJonesParameters(1.0, 1.0, 2.5)

    system = PotentialNBodySystem(bodies, Dict(:lennard_jones => potential_parameters))
    boundary_conditions = CubicPeriodicBoundaryConditions(L)
    sim = NBodySimulation(system, (0.0, steps * Δt), boundary_conditions, NBodySimulator.NullThermostat(), 1.0)
    @time result = run_simulation(sim, integrator(), dt = Δt)

    return result
end

function sim_lennard_jones_fluid(sim::NBodySimulation,Δt,n_steps,integrator)
    sim_alt=NBodySimulation(sim.system, (0.0, n_steps * Δt), sim.boundary_conditions, sim.thermostat,sim.kb)
    @time result= run_simulation(sim_alt, integrator(),dt=Δt)
    return result
end