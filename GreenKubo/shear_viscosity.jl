using Molly, LinearAlgebra

T=1.0
ρ=0.7
dt=1e-3
γ=1.0
t_eq=100.0
t_corr=2.0
Npd=10
splitting="BAOAB"
r_c=2.5
output_file="autocorrelations_sv.out"

include("../utils/place_atoms.jl")

N=Npd^3
L=(N/ρ)^(1//3)
boundary=CubicBoundary(L,L,L)
@assert 2r_c<= L "Cutoff radius too large relative to domain"

nf = CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff= 1.2r_c,unit_cell=boundary)


atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,boundary)
velocities=[velocity(1.0,T,1.0) for i=1:N]
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))
l_corr=ceil(Int64,t_corr/dt)
n_iter_sim=20000

sim=LangevinSplitting(dt=dt,temperature=T,friction=γ,splitting="BAOAB",remove_CM_motion=false)
norm_cst=inv(sqrt(3N))

function R(sys,args...;kwargs ...)
    N=length(sys)
    x_velocities=view(reinterpret(reshape,Float64,sys.velocities),1,:)
    y_coords=view(reinterpret(reshape,Float64,sys.coords),2,:)
    return dot(x_velocities,sin.(2π*y_coords/L))
end

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),neighbor_finder=nf,boundary=boundary,energy_units=NoUnits,force_units=NoUnits,k=1.0)
simulate!(sys,sim,n_steps_eq)
sys=System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),neighbor_finder=nf,boundary=boundary,energy_units=NoUnits,force_units=NoUnits,k=1.0,loggers=(autocorrelation=AutoCorrelationLogger(R,Float64,1,l_corr),))

for i=1:n_iter_sim
simulate!(sys,sim,n_steps_eq)
    f=open(output_file,"a")
    println(f,values(sys.loggers.autocorrelation;normalize=false))
    println(f,sys.loggers.autocorrelation.n_timesteps)
    close(f)
end

