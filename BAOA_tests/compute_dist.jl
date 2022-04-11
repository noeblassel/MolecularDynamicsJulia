using Random, LinearAlgebra

include("integrators.jl")
include("hist_utils.jl")
include("potential.jl")


println("usage: dt N_particles N_integration_steps N_equilibration_steps log_frequency splitting N_bins potential=PERIODIC|QUADRATIC|DOUBLE_WELL")

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential)
force_dict=Dict("PERIODIC"=>minus_d_periodic_potential,"QUADRATIC"=>minus_d_quadratic_potential,"DOUBLE_WELL"=>minus_d_double_well_potential)

dt=parse(Float64,ARGS[1])
N=parse(Int64,ARGS[2])
N_steps=parse(Int64,ARGS[3])
N_eq=parse(Int64,ARGS[4])
log_every=parse(Int64,ARGS[5])
splitting=ARGS[6]
N_bins=parse(Int64,ARGS[7])

potential=potential_dict[ARGS[8]]
force=potential_dict[ARGS[9]]

ps=zeros(N)
qs=zeros(N)


L=1.0
qlims = (ARGS[7]=="PERIODIC" ) ? (0.0,L) : (-10.0,10.0) 
plims= (-10.0,10.0)


bc = (ARGS[7]=="PERIODIC" ) ? PeriodicBoundaryCondition(L) : InfiniteBox()
sim=LangevinSplitting(dt=dt,γ=1.0,T=1.0,splitting=splitting,bc=bc)

hist=zeros(Int64,N_bins,N_bins)

simulate!(ps,qs,potential,force,hist,qlims,plims,sim,N_eq)
hist=zeros(Int64,N_bins,N_bins)


for i=1:log_every:N_steps
    simulate!(ps,qs,potential,force,hist,qlims,plims,sim,log_every)
    write_hist2D(hist,"bins_$(splitting)_$(ARGS[7])_$(dt).out")
end