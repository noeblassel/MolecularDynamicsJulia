using Statistics

include("../molly/MollyExtend.jl")


(ρmin,dρ,ρmax)=parse.(Float64,ARGS)
ρs=ρmin:dρ:ρmax

T = 1.2848906454490823

file=open("pressures.txt","a")

Npd = 12
N = Npd^3

r_c=4.0

atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

for ρ in ρs

    L = (N / ρ)^(1 // 3)

    box_size = SVector(L, L, L)

    coords = place_atoms_on_lattice(Npd, box_size)

    velocities = [reduced_velocity_lj(T) for i in 1:N]

    inter = LennardJones(cutoff = DistanceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

    nf = nothing
    nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)



    γ = 1.0
    eq_steps = 10000
    samp_steps = 30000

    dt = 5e-3

    seed = 1234

    simulator = VelocityVerlet(dt = dt)

    sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)

    simulate!(sys, simulator, eq_steps)

    loggers = Dict(:pressure => PressureLoggerReduced(Float64, 1))

    sys = System(atoms = sys.atoms, coords = sys.coords, velocities = sys.velocities, pairwise_inters = (inter,), box_size = sys.box_size, neighbor_finder = sys.neighbor_finder, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

    simulate!(sys, simulator, samp_steps)
    println(file,ρ," ", mean(sys.loggers[:pressure].pressures))
end

close(file)