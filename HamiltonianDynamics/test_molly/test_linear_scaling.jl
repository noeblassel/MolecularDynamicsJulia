using Molly,ProgressMeter,Plots

include("../molly/SimulateLennardJones.jl")
include("../molly/custom_simulators.jl")
include("../molly/custom_loggers.jl")

times=[]
ρ=0.5
T=0.5

for Npd in 1:12
    t=@elapsed sim_lennard_jones_fluid(Npd,ρ ,T,5e-3,5000,VelocityVerlet,[],3.0)
    push!(times,t)

    
end
plot((1:12).^3,times,xlabel="N",ylabel="runtime",label="",title="Verlet,no neighbor finder")