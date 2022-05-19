include("../molly/MollyExtend.jl")
using .MollyExtend

println("Usage: T ρ dt γ η forcing_type=COLOR|SINGLE t_equilibration t_simulation N_atoms_per_dimension scheme cutoff_radius")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
η=parse(Float64,ARGS[5])
forcing_type=ARGS[6]
t_eq=parse(Float64,ARGS[7])
t_sim=parse(Float64,ARGS[8])

Npd=parse(Int64,ARGS[9])
splitting=ARGS[10]
r_c=parse(Float64,ARGS[11])

N=Npd^3
L=(N/ρ)^(1//3)
box_size=SVector(L,L,L)

@assert 2r_c<= L "Cutoff radius too large relative to domain"

nf = (3.6r_c < L) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=1.2r_c)


atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=init_velocities(T,[a.mass for a=atoms],1.0)
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)
(forcing_type == "COLOR") && (η*=inv(N))

forcing = (forcing_type== "COLOR") ? ColorDriftNEMD(N,η,3) : SingleDriftNEMD(N,1,η)


ff=forcing.force_field
R= (forcing_type== "COLOR") ? MobilityObservable(ff) : (s::System,neighbors=nothing) -> s.velocities[1][1]

n_steps_eq=Int64(floor(t_eq/dt))
n_steps_sim=Int64(floor(t_sim/dt))

sim_eq=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
sim=LangevinSplitting(dt=dt,γ=γ,T=T,splitting=splitting)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),general_inters=(forcing,),box_size=box_size,neighbor_finder=nf,loggers=Dict{Symbol,Any}(),force_units=NoUnits,energy_units=NoUnits)

simulate!(sys,sim_eq,n_steps_eq)
sys.loggers=Dict(:mobility=>GeneralObservableLogger(R,1),:elapsed_time=>ElapsedTimeLogger(),:meta=>LogLogger([:mobility,:elapsed_time],["mobility_estimates$(η)_$(forcing_type).out","elapsed_time.out"],[100_000,10000],[false,true],["a","a"]))
simulate!(sys,sim,n_steps_sim)