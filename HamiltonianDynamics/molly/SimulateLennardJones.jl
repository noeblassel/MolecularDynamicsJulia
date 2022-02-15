
include("custom_loggers.jl")
include("../utils/PlaceAtoms.jl")

log_dict = Dict(
    :position => n -> CoordinateLogger(Float64, n),
    :temperature => n -> TemperatureLogger(Float64, n),
    :hamiltonian => n -> HamiltonianLogger(Float64, n),
    :kinetic_energy => n -> KineticEnergyLoggerNoDims(Float64, n),
    :potential_energy => n -> PotentialEnergyLogger(Float64, n),
    :velocity => n -> VelocityLogger(Float64, n),
    :pressure => n -> PressureLoggerReduced(Float64, n),
    :virial => n -> VirialLogger(Float64, n)
)


"""
Runs a simulation for the Lennard-Jones fluid of the given species, with N_per_dim^3 particles at density ρ and temperature T.
If T and ρ are of a Real subtype and reduced_units=false, they must r   espectively be given in K and in nm^-3 (number of particles per nanometer).
Alternatively, they can be specified in physical (Unitful) units,or in reduced units with reduced_units set to true.
The parameter Δt is always dimensionless.
"""

function sim_lennard_jones_fluid(N_per_dim::Integer, ρ::Real, T::Real, Δt::Real, steps::Integer, integrator, observables, r_c::Real; equilibration_steps = 0)

    N = N_per_dim^3
    atoms = [Atom(mass = 1.0, σ = 1.0, ϵ = 1.0) for i in 1:N]#in reduced units
    L = (N / ρ)^(1 // 3)

    domain = SVector(L, L, L)

    coords = place_atoms_on_lattice(N_per_dim, domain)
    velocities = [reduced_velocity_lj(T) for i in 1:N]
    interactions = (LennardJones(cutoff = ShiftedPotentialCutoff(r_c), force_units = NoUnits, energy_units = NoUnits),)

    loggers = Dict()
    for (ob, n...) = observables
        if ob == :state
            loggers[ob] = StateLogger(n...)
        else
            loggers[ob] = log_dict[ob](first(n))
        end
    end

    if equilibration_steps>0
        simulator = integrator(dt = Δt, coupling = RescaleThermostat(T))
        sys_eq = System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, energy_units = NoUnits, force_units = NoUnits)
        simulate!(sys_eq, simulator, equilibration_steps)
        sys = System(atoms = atoms, general_inters = interactions, coords = sys_eq.coords, velocities = sys_eq.velocities, box_size = domain, loggers = loggers, energy_units = NoUnits, force_units = NoUnits)
    else
        sys=System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, energy_units = NoUnits, force_units = NoUnits)
    end

    simulator = integrator(dt = Δt)
    @time simulate!(sys, simulator, steps)
    return sys
end


function sim_lennard_jones_fluid!(sys::System,Δt::Real,steps::Integer,integrator)
    simulator=integrator(dt=Δt)
    @time simulate!(sys,simulator,steps)
    return sys
end

reduced_velocity_lj(T::Real) = sqrt(T::Real) * (@SVector randn(3))
