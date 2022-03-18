include("../molly/MollyExtend.jl")

seed = 2022
Npd = 6
N = Npd^3
ρ = 0.4
L = (N / ρ)^(1 // 3)
T = 1.5
dt = 5e-3
r_c = 4.0

n_steps = 5000

box_size = SVector(L, L, L)
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)

nf = DistanceNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c)
coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T, atoms[i].mass) for i in 1:N]

sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), box_size=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=Dict{Symbol,Any}())
γ=1.0
sim= LangevinSplitting(dt=dt, γ=γ, T=T, splitting="BOBOBABABAOBAB", rseed=seed)
simulate!(sys,sim,100)
@time simulate!(sys,sim,5000)