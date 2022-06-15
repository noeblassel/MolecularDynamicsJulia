include("../molly/MollyExtend.jl")

using .MollyExtend

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3
γ=1.0
dt=5e-3
v=1.0

L = (N / ρ)^(1 // 3)
box_size = SVector(L, L, L)

coords = place_atoms_on_3D_lattice(Npd,box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = init_velocities(T,[a.mass for a=atoms],1.0)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

G(y::Float64)=sin(2π*y/L)

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
simulator=NortonShearViscosityTest(dt = dt, γ = γ, T = T,v=v,G=G)

n_eq_steps=10
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,),box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)
#@profview simulate!(sys,simulator,n_eq_steps)