using Plots
using Molly
using ProgressMeter
using Statistics

include("../molly/input.jl")
include("../molly/KineticEnergyLoggerNoDims.jl")
include("../molly/SymplecticEuler.jl")

initial_system=read_reduced_lj_state("HamiltonianDynamics/test_molly/equilibrium.out")

loggers=Dict(:hamiltonian=>TotalEnergyLogger(Float64,1))
initial_system.loggers=loggers

stds=[]

N=length(initial_system)

energy_plot=plot(xlabel="time step",ylabel="H(t)/N")
dt_range=[0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064]

for dt in dt_range
    sys=deepcopy(initial_system)
    simulate!(sys,VelocityVerlet(dt=dt),1000)
    plot!(energy_plot,sys.loggers[:hamiltonian].energies/N,label="Δt=$(dt)")
    push!(stds,std(sys.loggers[:hamiltonian].energies))
end

savefig(energy_plot,"energy_plot.png")
savefig(plot(dt_range,stds,label="std of hamiltonian",title="Verlet energy conservation",xlabel="Δt",ylabel="std"),"stds.png")

