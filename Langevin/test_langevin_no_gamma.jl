using Plots, Statistics,LinearAlgebra
using Test
include("../molly/MollyExtend.jl")

ρ = 0.2
T = 1.0

r_c = 4.0

Npd = 10
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]

n_steps=2000

dt = 5e-3

seed = 1234

sim_langevin = LangevinTest(dt = dt, γ = 1.0, T = T, rseed = seed)
sim_verlet=VelocityVerlet(dt=dt)

r_c=4.0

nf_langevin=NoNeighborFinder()
nf_verlet=NoNeighborFinder()

lj=LennardJones(cutoff=DistanceCutoff(r_c),force_units=NoUnits,energy_units=NoUnits)

sys_langevin = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (lj,), box_size = box_size,neighbor_finder=nf_langevin, force_units = NoUnits, energy_units = NoUnits)
sys_verlet = System(atoms = deepcopy(atoms), coords = deepcopy(coords), velocities = deepcopy(velocities), pairwise_inters = (lj,), box_size = box_size,neighbor_finder=nf_verlet, force_units = NoUnits, energy_units = NoUnits)

@time simulate!(sys_langevin,sim_langevin,n_steps)
@time simulate!(sys_verlet,sim_verlet,n_steps)

"""for i=1:N
    @test sys_langevin.coords[i]==sys_verlet.coords[i]
end"""