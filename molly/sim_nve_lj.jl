log_dict = Dict(
    :position => n -> CoordinateLogger(Float64, n),
    :temperature => n -> TemperatureLoggerReduced(Float64, n),
    :hamiltonian => n -> HamiltonianLogger(Float64, n),
    :kinetic_energy => n -> KineticEnergyLoggerNoDims(Float64, n),
    :potential_energy => n -> PotentialEnergyLogger(Float64, n),
    :velocity => n -> VelocityLogger(Float64, n),
    :pressure => n -> PressureLoggerReduced(Float64, n),
    :virial => n -> VirialLogger(Float64, n),
    :state => (n,fp)->StateLogger(n,fp)
)


"""
Runs a simulation for the Lennard-Jones fluid of the given species, with N_per_dim^3 particles at density ρ and temperature T.
If T and ρ are of a Real subtype and reduced_units=false, they must r   espectively be given in K and in nm^-3 (number of particles per nanometer).
Alternatively, they can be specified in physical (Unitful) units,or in reduced units with reduced_units set to true.
The parameter Δt is always dimensionless.
"""

function sim_lennard_jones_fluid_nve(N_per_dim::Integer, ρ::Real, T_ini::Real, Δt::Real, steps::Integer, integrator, observables, r_c::Real; equilibration_steps = 0)

    N = N_per_dim^3
    atoms = [Atom(mass = 1.0, σ = 1.0, ϵ = 1.0) for i in 1:N]#in reduced units
    L = (N / ρ)^(1 // 3)

    domain = SVector(L, L, L)

    coords = place_atoms_on_lattice(N_per_dim, domain)
    velocities = [reduced_velocity_lj(T) for i in 1:N]
    interactions = (LennardJones(cutoff = DistanceCutoff(r_c), force_units = NoUnits, energy_units = NoUnits, nl_only = true),)
    nb_mat = trues(N, N)
    mat_14 = falses(N, N)

    nf = nothing

    if L<3r_c
        nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)##in case the system is too small, revert to tree neighbor finder
    else
        nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
    end


    loggers = Dict()
    
    for (ob, n...) = observables
        loggers[ob]=log_dict[ob](n...)
    end

    if equilibration_steps > 0
        simulator = integrator(dt = Δt)
        sys_eq = System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf)
        simulate!(sys_eq, simulator, equilibration_steps)
        sys = System(atoms = atoms, general_inters = interactions, coords = sys_eq.coords, velocities = sys_eq.velocities, box_size = domain, loggers = loggers, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf)
    else
        sys = System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, loggers = loggers, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf)
    end

    simulator = integrator(dt = Δt)
    @time simulate!(sys, simulator, steps)
    return sys
end


function sim_lennard_jones_fluid_nve!(sys::System, Δt::Real, steps::Integer, integrator)
    simulator = integrator(dt = Δt)
    @time simulate!(sys, simulator, steps)
    return sys
end