using NBodySimulator,StaticArrays,Plots,Statistics

include("../nbodysimulator/io.jl")
include("../nbodysimulator/custom_simulators.jl")
include("../nbodysimulator/sim_nve_lj.jl")
include("../../utils/reduced_units.jl")

N_pd=6

ρ=0.5
T=0.5

res=sim_lennard_jones_fluid(N_pd,ρ,T,0.005,5000,VelocityVerlet)
save_reduced_lj_state(res,"ini_config.out")
sim=read_reduced_lj_state("ini_config.out")
stds=[]

N=length(res.simulation.system.bodies)
energy_plot=plot(xlabel="time step",ylabel="H(t)/N")
dt_range=[0.0001,0.0004,0.0008,0.0016,0.0032,0.005,0.0064]

for dt in dt_range
    res=sim_lennard_jones_fluid(sim,dt,1000,VelocityVerlet)
    energies=[total_energy(res,t) for t in res.solution.t]
    plot!(energy_plot,energies/N,label="Δt=$(dt)")
    push!(stds,std(energies))
end

savefig(energy_plot,"energy_plot.png")
savefig(plot(dt_range,stds,label="std of hamiltonian",title="Verlet energy conservation",xlabel="Δt",ylabel="std"),"stds.png")

