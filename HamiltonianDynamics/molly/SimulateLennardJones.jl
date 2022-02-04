using Molly

include("SymplecticEuler.jl")
include("PressureLogger.jl")
include("HamiltonianLogger.jl")
include("KineticEnergyLoggerNoDims.jl")
include("../ReducedUnits.jl")
include("PlaceAtoms.jl")

const log_dict = Dict(:position => n->CoordinateLogger(Float64,n), :pressure => n->PressureLogger(Float64,n), :temperature => n->TemperatureLogger(Float64,n), :hamiltonian => n->HamiltonianLogger(Float64,n),:kinetic_energy=>n->KineticEnergyLoggerNoDims(Float64,n),:potential_energy=>n->PotentialEnergyLogger(Float64,n),:velocity=>n->VelocityLogger(Float64,n))


"""
Runs a simulation for the Lennard-Jones fluid of the given species, with N_per_dim^3 particles at density ρ and temperature T.
If T and ρ are of a Real subtype and reduced_units=false, they must respectively be given in K and in nm^-3 (number of particles per nanometer).
Alternatively, they can be specified in physical (Unitful) units,or in reduced units with reduced_units set to true.
The parameter Δt is always dimensionless.
"""

function sim_lennard_jones_fluid(N_per_dim, ρ, T, Δt, steps,integrator, observables;species=nothing,reduced_units=true)
    
    if reduced_units
        ρʳ=ρ
        Tʳ=T
    else 
        @assert species!==nothing "cannot compute reduced units for unspecified species"
        ρʳ=get_reduced_density(species,ρ)
        Tʳ=get_reduced_temperature(species,T)
    end

    N = N_per_dim^3
    atoms = [Atom(mass = 1.0, σ = 1.0, ϵ = 1.0) for i in 1:N]
    L = (N/ ρʳ)^(1 // 3)

    domain = SVector(L, L, L)

    coords=place_atoms_on_lattice(N_per_dim,domain)
    velocities = [velocity(1.0, Tʳ) for i in 1:N]
    interactions = (LennardJones(cutoff = DistanceCutoff(2.5),force_units=NoUnits,energy_units=NoUnits),)

    loggers = Dict()
    for (ob, n) = observables
        loggers[ob] = log_dict[ob](n)
    end
    sys = System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, loggers = loggers,energy_units=NoUnits,force_units=NoUnits)
    simulator = integrator(dt = Δt)
    simulate!(sys, simulator, steps)
    return sys
end
