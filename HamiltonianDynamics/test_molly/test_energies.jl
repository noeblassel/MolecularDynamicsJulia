using Plots
using Statistics

include("../../molly/MollyExtend.jl")

initial_system=read_reduced_lj_state("starting_states/config.lj_out")

loggers=Dict(:hamiltonian=>TotalEnergyLogger(Float64,1))
initial_system.loggers=loggers

spreads=[]

N=length(initial_system)

#energy_plot=plot(xlabel="time step",ylabel="H(t)")
dt_range=0.0001:0.0001:0.0064

for dt in dt_range
    sys=deepcopy(initial_system)
    simulate!(sys,SymplecticEulerB(dt=dt),500)
    #plot!(energy_plot,sys.loggers[:hamiltonian].energies,label="Δt=$(dt)")
    push!(spreads,maximum(sys.loggers[:hamiltonian].energies)-minimum(sys.loggers[:hamiltonian].energies))
end

f=open("energy_fluctuations_SEB.txt","w")

for (dt,dH) in zip(dt_range,spreads)
    println(f,dt," ",dH)
end

close(f)

#savefig(energy_plot,"energy_fluctuations_verlet.pdf")
#savefig(plot(dt_range,stds,label="std of hamiltonian",title="Verlet energy conservation",xlabel="Δt",ylabel="std"),"stds.png")

