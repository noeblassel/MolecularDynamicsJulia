using Random, LinearAlgebra

include("integrators.jl")
include("hist_utils.jl")
include("potential.jl")


println("usage: dt N_particles N_integration_steps N_equilibration_steps log_frequency splitting N_bins potential=PERIODIC|QUADRATIC|DOUBLE_WELL|TILTED_DOUBLE_WELL γ")
force_dict=Dict("PERIODIC"=>grad_periodic_potential,"QUADRATIC"=>grad_quadratic_potential,"DOUBLE_WELL"=>grad_double_well_potential,"TILTED_DOUBLE_WELL"=>grad_tilted_double_well_potential)

dt=parse(Float64,ARGS[1])
N=parse(Int64,ARGS[2])
N_steps=parse(Int64,ARGS[3])
N_eq=parse(Int64,ARGS[4])
log_every=parse(Int64,ARGS[5])
splitting=ARGS[6]
N_bins=parse(Int64,ARGS[7])
∇V=force_dict[ARGS[8]]
γ=parse(Float64,ARGS[9])

∇K(p::Float64)=sign(p)*p^2

ps=zeros(N)
qs=zeros(N)

L=1.0
qlims = (ARGS[8]=="PERIODIC" ) ? (0.0,L) : (-5.0,5.0) 
plims= (-5.0,5.0)


bc = (ARGS[8]=="PERIODIC" ) ? PeriodicBoundaryCondition(L) : InfiniteBox()
sim=LangevinSplitting(dt=dt,γ=γ,T=1.0,splitting=splitting,bc=bc)

hist=zeros(Int64,N_bins,N_bins)

simulate!(ps,qs,∇V,∇K,hist,qlims,plims,sim,N_eq)
hist=zeros(Int64,N_bins,N_bins)

for i=1:log_every:N_steps
    simulate!(ps,qs,∇V,∇K,hist,qlims,plims,sim,log_every)
    write_hist2D(hist,"bins_$(splitting)_$(ARGS[8])_$(dt)_$(γ).out")
end
