using Molly,LinearAlgebra,Random

include("../molly/MollyExtend.jl")
include("integrators.jl")

using .MollyExtend

metropolis=parse(Bool,ARGS[1])

Npd=4
N=Npd^3

dt=5e-7


ρ=0.2
L = (N / ρ)^(1 // 3)
box_size=SVector(L,L,L)

atoms=[Atom(ϵ=1.0,mass=1.0,σ=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=init_velocities(1.0,[a.mass for a=atoms],1.0)

#------------------------   ----------------------------------
r_c_nf=2.5

if (L > 3 * r_c_nf) && N>900
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf, unit_cell = box_size)
elseif N > 900
    nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf)
else
    nf = DistanceNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf)
end
#-----------------------   -----------------------------------
r_c_inter=2.0

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c_inter), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = Dict{Symbol,Any}())

n_steps_eq=100_000
sim=MALA(dt=1e-6,T=1.0)
simulate!(sys,sim,n_steps_eq)


n_steps_sim=10_000_000
for lg_dt in range(-7,-5,20)
    dt=10^lg_dt
    sim=MALA(dt=dt,T=1.0,is_metropolis=metropolis)
    println(dt," ",sim.n_accepted/sim.n_total)
    flush(stdout)
end