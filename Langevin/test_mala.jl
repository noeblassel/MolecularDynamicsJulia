using Molly

include("../molly/MollyExtend.jl")

using .MollyExtend

Npd = 3
N = Npd^3
ρ = 0.4

T = 1.0
L = (N / ρ)^(1 // 3)
lg_dts = range(-5, -3, 40)
Tfin = parse(Float64, ARGS[1])
proposal = ARGS[2]
is_metropolis = (ARGS[3] == "METROPOLIS")

dt_eq = 1e-3
N_steps_eq = 1_000_000
dts = 10 .^ lg_dts

atoms = [Atom(ϵ=1.0, σ=1.0, mass=1.0) for i = 1:N]
box_size = SVector(L, L, L)
coords = place_atoms_on_3D_lattice(Npd, box_size)
velocities = init_velocities(T, mass.(atoms), 1.0)

inters = (LennardJones(energy_units=NoUnits, force_units=NoUnits),)
sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=inters, loggers=Dict{Symbol,Any}(), box_size=box_size, force_units=NoUnits, energy_units=NoUnits)

obs1(s::System, neighbors=Nothing) = potential_energy(s, neighbors) / N
obs2(s::System, neighbors=Nothing) = (X = forces(s, neighbors); dot(X, X) / N)

Vs = Float64[]
grad_Vs = Float64[]
Rs = Float64[]
println("equilibriating")
eq_sim = MALA(dt=dt_eq, T=T)
simulate!(sys, eq_sim, N_steps_eq)
f = open("mala_output_$(proposal)_$(ARGS[3]).out", "w")
println(f, "dt V grad_V R")

for dt in reverse(dts)
    println(dt)
    sys.loggers = Dict(:potential => AverageObservableLogger(obs1, 1), :grad => AverageObservableLogger(obs2, 1))

    N_samples = round(Int64, Tfin / dt)

    if proposal == "EM"
        sim = MALA(dt=dt, T=T,is_metropolis=is_metropolis)
    else
        sim = MALA_HMC(dt=dt, T=T,is_metropolis=is_metropolis)
    end

    simulate!(sys, sim, N_samples)
    println(f, "$(dt) $(sys.loggers[:potential].sum/sys.loggers[:potential].n_samples) $(sys.loggers[:grad].sum/sys.loggers[:potential].n_samples) $(1-sim.n_accepted/sim.n_total)")
    flush(f)

end
close(f)