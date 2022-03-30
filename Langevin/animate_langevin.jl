using Plots,LinearAlgebra
include("../molly/MollyExtend.jl")
include("../utils/animate.jl")
ρ = 0.25

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
T = [(i<= N/2) ? 8.0 : 0 for i=1:N]
velocities = [reduced_velocity_lj(T[i],atoms[i].mass) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing

try
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
catch
    println("revert to tree nf")
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

loggers = Dict(:coords => CoordinateLogger(Float64, 1))

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)


γ = 1.0
n_steps = 20000
dt = 1e-2

seed = 1234

simulator = LangevinSplitting(dt = dt, γ = γ, T = T, rseed = seed,splitting="BAOAB")
@time simulate!(sys, simulator, n_steps)
animate_trajectories(sys.loggers[:coords].coords,"langevin_inhomogenous_T.mp4";ix=N,camera=(0,0))