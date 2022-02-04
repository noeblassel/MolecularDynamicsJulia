using Molly

include("SymplecticEuler.jl")
include("PressureLogger.jl")
include("ReducedUnits.jl")


const r_c=3σ
const log_dict = Dict("positions" => n->CoordinateLogger(Float64,n), "pressure" => n->PressureLogger(Float64,n), "temperature" => n->TemperatureLogger(Float64,n), "hamiltonian" => n->TotalEnergyLogger(Float64,n),"kinetic_energy"=>n->KineticEnergyLogger(Float64,n),"potential_energy"=>n->PotentialEnergyLogger(Float64,n))


function sim_argon(N_per_dim, ρ, T, Δt, sim_duration,integrator, observables)
    steps = Int32(ceil(sim_duration / Δt))
    N = N_per_dim^3
    atoms = [Atom(mass = m, σ = σ, ϵ = ϵ) for i in 1:N]
    L = (N * m / ρ)^(1 // 3)
    l = L / N_per_dim
    domain = SVector(L, L, L)

    coords = [SVector(i * l, j * l, k * l) for i = 1:N_per_dim, j = 1:N_per_dim, k = 1:N_per_dim]#arrange atoms on an orhogonal lattice
    coords = reshape(coords, N)

    velocities = [velocity(m, T) for i in 1:N]
    interactions = (LennardJones(cutoff = DistanceCutoff(r_c),force_units=NoUnits,energy_units=NoUnits),)

    loggers = Dict()
    for (ob, n) = observables
        loggers[ob] = log_dict[ob](n)
    end
    sys = System(atoms = atoms, general_inters = interactions, coords = coords, velocities = velocities, box_size = domain, loggers = loggers,energy_units=NoUnits,force_units=NoUnits)
    simulator = integrator(dt = Δt)
    simulate!(sys, simulator, steps)
    return sys
end
