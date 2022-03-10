include("../molly/MollyExtend.jl")

println("usage: generate_state.jl T ρ dt r_c Npd n_steps output_file")
T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
r_c=parse(Float64,ARGS[4])
Npd=parse(Int64,ARGS[5])
n_steps=parse(Int64,ARGS[6])
output_file=ARGS[7]


N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf=nothing

if L>3*r_c
    nf=CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c,unit_cell=box_size)
else
    nf=TreeNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c)
end
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)


γ=1.0
simulator=LangevinBAOAB(T=T,γ=γ,dt=dt)
simulate!(sys,simulator,n_steps)

save_reduced_lj_state(sys,output_file)