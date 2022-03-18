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
sys_abstract = deepcopy(sys)

γ = 1.0

sim_bao = LangevinBAO(dt=dt, γ=γ, T=T, rseed=seed)
sim_bao_abstract = LangevinSplitting(dt=dt, γ=γ, T=T, splitting="BAO", rseed=seed)

sim_babo = LangevinBABO(dt=dt, γ=γ, T=T, rseed=seed)
sim_babo_abstract = LangevinSplitting(dt=dt, γ=γ, T=T, splitting="BABO", rseed=seed)

sim_baoa = LangevinBAOA(dt=dt, γ=γ, T=T, rseed=seed)
sim_baoa_abstract = LangevinSplitting(dt=dt, γ=γ, T=T, splitting="BAOA", rseed=seed)

sim_baoab = LangevinBAOAB(dt=dt, γ=γ, T=T, rseed=seed)
sim_baoab_abstract = LangevinSplitting(dt=dt, γ=γ, T=T, splitting="BAOAB", rseed=seed)

println("-----BAO-----")

simulate!(sys, sim_bao, n_steps)
simulate!(sys_abstract, sim_bao_abstract, n_steps)


@time simulate!(sys, sim_bao, n_steps)
@time simulate!(sys_abstract, sim_bao_abstract, n_steps)

@assert all(sys.velocities[i] ≈ sys_abstract.velocities[i] for i = 1:N)
@assert all(sys.coords[i] ≈ sys_abstract.coords[i] for i = 1:N)

println("-----BAOA-----")

simulate!(sys, sim_baoa, n_steps)
simulate!(sys_abstract, sim_baoa_abstract, n_steps)


@time simulate!(sys, sim_baoa, n_steps)
@time simulate!(sys_abstract, sim_baoa_abstract, n_steps)

@assert all(sys.velocities[i] ≈ sys_abstract.velocities[i] for i = 1:N)
@assert all(sys.coords[i] ≈ sys_abstract.coords[i] for i = 1:N)

println("-----BABO-----")

simulate!(sys, sim_babo, n_steps)
simulate!(sys_abstract, sim_babo_abstract, n_steps)


@time simulate!(sys, sim_babo, n_steps)
@time simulate!(sys_abstract, sim_babo_abstract, n_steps)

@assert all(sys.velocities[i] ≈ sys_abstract.velocities[i] for i = 1:N)
@assert all(sys.coords[i] ≈ sys_abstract.coords[i] for i = 1:N)

println("-----BAOAB-----")

simulate!(sys, sim_baoab, n_steps)
simulate!(sys_abstract, sim_baoab_abstract, n_steps)


@time simulate!(sys, sim_baoab, n_steps)
@time simulate!(sys_abstract, sim_baoab_abstract, n_steps)

@assert all(sys.velocities[i] ≈ sys_abstract.velocities[i] for i = 1:N)
@assert all(sys.coords[i] ≈ sys_abstract.coords[i] for i = 1:N)