using Pkg, Statistics

Pkg.instantiate()

include("../molly/MollyExtend.jl")



#julia test_bias.jl T ρ dt tfin simulator output

@assert length(ARGS) == 10 "Error (Wrong Argument Count) Usage: test_bias.jl T ρ Δt teq tfin  Npd Nruns r_c splitting OUTPUT|STDOUT"
T = parse(Float64, ARGS[1])
ρ = parse(Float64, ARGS[2])
dt = parse(Float64, ARGS[3])
teq = parse(Float64, ARGS[4])
tfin = parse(Float64, ARGS[5])
Npd = parse(Int64, ARGS[6])
Nruns = parse(Int64, ARGS[7])
r_c = parse(Float64, ARGS[8])
sim = ARGS[9]
output_file = ARGS[10]

if !isfile(output_file)
    g = output_file!="STDOUT" ? open(output_file, "w") : stdout
    println(g, "Trajectorial averages of NVT Lennard-Jones system of $(Npd^3) particles at T=$(T), ρ=$(ρ) with $(sim) splitting, $(r_c) shifted force cutoff. Physical time of each run: $(tfin). All units are reduced.")
    println(g, "[dt] [average potential energy] [average kinetic energy] [average virial]")
    (output_file!="STDOUT") && close(g)
end

dt_eq = 5e-3
eq_nsteps = Int64(round(teq / dt_eq))
N = Npd^3
L = (N / ρ)^(1 // 3)
box_size = SVector(L, L, L)
inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing
n_steps = Int64(round(tfin / dt))

if (L > 3 * r_c) && N>900
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
elseif N > 900
    nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
else
    nf = DistanceNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T, atoms[i].mass) for i in 1:N]
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = Dict{Symbol,Any}())

γ = 1.0

simulator = LangevinSplitting(T = T, γ = γ, dt = dt_eq,splitting=sim)
simulate!(sys, simulator, eq_nsteps)
loggers = Dict(:potential_energy => PotentialEnergyLogger(Float64, 1), :kinetic_energy => KineticEnergyLoggerNoDims(Float64, 1), :virial => VirialLogger(Float64, 1))
sys.loggers = loggers

for i = 1:Nruns

    simulate!(sys, simulator, n_steps)

    f = output_file !="STDOUT" ? open(output_file, "a") : stdout
    Vhat = mean(sys.loggers[:potential_energy].energies)
    Khat = mean(sys.loggers[:kinetic_energy].energies)
    What = mean(sys.loggers[:virial].energies)

    println(f, "$(dt) $(Vhat) $(Khat) $(What)")
    (output_file!="STDOUT") && close(f)

    empty!(sys.loggers[:potential_energy].energies)
    empty!(sys.loggers[:kinetic_energy].energies)
    empty!(sys.loggers[:virial].energies)
end