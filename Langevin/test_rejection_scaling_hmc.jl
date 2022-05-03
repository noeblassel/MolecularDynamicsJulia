include("../molly/MollyExtend.jl")
using .MollyExtend

@assert length(ARGS) == 5 "Error(Wrong Argument Count) expected 5, got $(length(ARGS)).\n Usage: test_rejection_scaling_hmc.jl T ρ teq tfin Npd"
T = parse(Float64, ARGS[1])
ρ = parse(Float64, ARGS[2])
teq = parse(Float64, ARGS[3])
tfin = parse(Float64, ARGS[4])
Npd = parse(Int64, ARGS[5])
r_c_inter=1.6


dt_eq = 5e-3
eq_nsteps = Int64(round(teq / dt_eq))
N = Npd^3
L = (N / ρ)^(1 // 3)
box_size = SVector(L, L, L)
inter = LennardJones(cutoff = ShiftedForceCutoff(r_c_inter), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

#nf = nothing
r_c=1.8
if (L > 3 * r_c) && N>900
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
elseif N > 900
    nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
else
    nf = DistanceNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

coords = place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity(T, atoms[i].mass,1.0) for i in 1:N]
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = Dict{Symbol,Any}())

γ = 1.0

simulator_eq= LangevinSplitting(T=T, γ=γ, dt= dt_eq, splitting="BAOAB")
simulate!(sys,simulator_eq,eq_nsteps)

lg_dt_range=range(-5,-4,20)
dts= 10 .^lg_dt_range
for dt in reverse(dts)
    n_steps = round(Int64,tfin / dt)    
    simulator = LangevinGHMC(T = T, γ = γ, dt = dt)
    simulate!(sys,simulator,n_steps)
    println(dt," ",1-simulator.n_accepted/simulator.n_total)
end