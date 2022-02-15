using NBodySimulator,BenchmarkTools,Plots

include("../utils/ReducedUnits.jl")
include("../nbodysimulator/SimulateLennardJones.jl")
include("../nbodysimulator/io.jl")

times=[]
T=0.5
ρ=0.5
steps=5000

for N=1:12
    t=@elapsed sim_lennard_jones_fluid(N,ρ,T,5e-3,1000,VelocityVerlet)
    push!(times,t)
    println(N)
end

plot([N^3 for N in 1:12],times,xlabel="N",ylabel="T",title="runtime scaling",label="")