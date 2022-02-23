using Molly, ProgressMeter

include("../molly/custom_simulators.jl")
include("../molly/custom_observables.jl")
include("../molly/custom_loggers.jl")
include("../molly/custom_cutoffs.jl")

include("../utils/ReducedUnits.jl")
include("../utils/PlaceAtoms.jl")
include("../utils/animate.jl")

ρ = 0.5
T = 0.5

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
velocities = [reduced_velocity_lj(T) for i in 1:N]
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff_fixed(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing

try
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
catch
    println("revert to tree nf")
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

loggers = Dict(:coords => CoordinateLogger(Float64, 1))

sys = System(atoms = atoms, coords = coords, velocities = velocities, general_inters = (), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)


γ = 1.0
n_steps = 20000
dt = 5e-3

seed = 1234

simulator = LangevinTest(dt = dt, γ = γ, T = T, rseed = seed)
@time simulate!(sys, simulator, n_steps)

using Plots,LinearAlgebra

animate_trajectories(sys.loggers[:coords].coords,"langevin_T0.1.mp4")