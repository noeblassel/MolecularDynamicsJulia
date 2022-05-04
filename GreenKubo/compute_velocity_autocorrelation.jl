T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
t_eq=parse(Float64,ARGS[5])
t_sim=parse(Float64,ARGS[6])
log_freq=parse(Float64,ARGS[7])

Npd=parse(Int64,ARGS[8])
splitting=ARGS[9]
r_c=parse(Float64,ARGS[10])
output_file=ARGS[11]

include("../molly/MollyExtend.jl")
using .MollyExtend,LinearAlgebra

println("Usage: T ρ dt γ t_equilibration t_simulation logging_frequency N_atoms_per_dimension scheme cutoff_radius")

N=Npd^3
L=(N/ρ)^(1//3)
box_size=SVector(L,L,L)

@assert 2r_c<= L "Cutoff radius too large relative to domain"

nf = (4.5r_c < L) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff= 1.5r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=1.5r_c)


atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=init_velocities(T,[a.mass for a=atoms],1.0)
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))
n_steps_sim=Int64(floor(t_sim/dt))
n_log_freq=Int64(floor(log_freq/dt))

sim_eq=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
sim=LangevinSplitting(dt=dt,γ=γ,T=T,splitting=splitting)

