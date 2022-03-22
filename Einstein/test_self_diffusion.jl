include("../molly/MollyExtend.jl")


T=2.5
ρ=0.7
Npd=5
N=Npd^3
L = (N / ρ)^(1 // 3)
r_c=4.0

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

loggers = Dict(:self_diffusion => SelfDiffusionLogger(1,coords;record_history=true))

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

sim=LangevinSplitting(dt=5e-3,γ=1.0,T=T,splitting="BAOAB")
simulate!(sys,sim,10000)