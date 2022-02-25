using Statistics

include("../molly/MollyExtend.jl")


ρmin=parse(Float64,ARGS[1])
ρs=ρmin:0.001:(0.0075125+ρmin)

T = 1.2848906454490823

file=open("pressures.txt","a")

Npd = 10
N = Npd^3

r_c=4.0

atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

for ρ in ρs

    L = (N / ρ)^(1 // 3)

    box_size = SVector(L, L, L)

    coords = place_atoms_on_lattice(Npd, box_size)

    velocities = [reduced_velocity_lj(T) for i in 1:N]

    inter = LennardJones(cutoff = DistanceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c,unit_cell=box_size)



    γ = 1.0
    eq_steps = 5000
    samp_steps = 10000

    dt = 5e-3

    simulator = LangevinTest(dt = dt,γ=γ,T=T)

    sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)

    simulate!(sys, simulator, eq_steps)

    loggers = Dict(:pressure => PressureLoggerNVT(T,Float64, 1))

    sys = System(atoms = sys.atoms, coords = sys.coords, velocities = sys.velocities, pairwise_inters = (inter,), box_size = sys.box_size, neighbor_finder = sys.neighbor_finder, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

    simulate!(sys, simulator, samp_steps)
    println(file,ρ," ", mean(sys.loggers[:pressure].pressures))
end

close(file)