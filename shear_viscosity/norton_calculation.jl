using Molly

include("../utils/place_atoms.jl")
include("norton_integrators.jl")
println("Usage: T ρ dt γ η forcing_type=SINUSOIDAL|LINEAR|CONSTANT t_equilibration n_iter_sim N_atoms_per_dimension cutoff_radius")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
η=parse(Float64,ARGS[5])
forcing_type=ARGS[6]
t_eq=parse(Float64,ARGS[7])
n_iter_sim=parse(Int64,ARGS[8])

Npd=parse(Int64,ARGS[9])
r_c=parse(Float64,ARGS[10])

N=Npd^3
L=(N/ρ)^(1//3)
box_size=CubicBoundary(L,L,L)

max_speed=10.0*sqrt(T)
n_steps_neighbors=floor(Int64,0.2*r_c/(dt*max_speed))


F_sin(y) = sin(2π * y / L)
F_const(y) = (y < L / 2) ? 1.0 : -1.0
F_lin(y) = (y < L / 2) ? 4 * (y - L / 4) / L : 4 * (3L / 4 - y) / L

G_imag(y) = sin(2π * y / L) / N
G_real(y) = cos(2π * y/ L) / N

F=F_sin
G=G_imag

if forcing_type=="LINEAR"
    F=F_lin
    G=G_real
elseif forcing_type=="CONSTANT"
    F=F_const
end

simulator=NortonSVIntegrator(dt,η,T,γ,F,G)

nf = (3.6r_c < L) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff=1.2r_c)

coords = place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=floor(Int64,t_eq/dt)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)

_=simulate!(sys,simulator,n_steps_eq)

for i=1:n_iter_sim
    force = simulate!(sys,simulator,n_steps_eq)
    f=open("norton_forcing_$(forcing_type)_$(η).out","a")
    write(f,force)
    close(f)
end