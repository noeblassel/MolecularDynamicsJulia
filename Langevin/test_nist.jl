
include("../molly/MollyExtend.jl")

ρs = 0.4502:0.005:0.5704

T = 1.2848906454490823


Npd = 12
N = Npd^3

for ρ in ρs

    L = (N / ρ)^(1 // 3)

    box_size = SVector(L, L, L)

    coords = place_atoms_on_lattice(Npd, box_size)

    velocities = [SVector(1.0, 1.0, 1.0) for i in 1:N]

    inter = LennardJones(cutoff = ShiftedForceCutoff_fixed(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

    nf = nothing

    try
        global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
    catch
        global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
    end


    γ = 1.0
    eq_steps = 10000
    samp_steps = 30000

    dt = 5e-3

    seed = 1234

    simulator = LangevinTest(dt = dt, γ = γ, T = T, rseed = seed)

    sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), general_inters = (), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)

    simulate!(sys, simulator, eq_steps)

    loggers = Dict(:pressure => PressureLoggerReduced(Float64, 1))

    sys = System(atoms = sys.atoms, coords = sys.coords, velocities = sys.velocities, general_inters = (), box_size = sys.box_size, neighbor_finder = sys.neighbor_finder, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

    simulate!(sys, simulator, samp_steps)
    println(mean(sys.loggers[:pressure].pressures))
end