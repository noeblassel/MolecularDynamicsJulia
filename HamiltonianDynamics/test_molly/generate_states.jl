using Molly
using ProgressMeter

include("../molly/io.jl")
include("../molly/SimulateLennardJones.jl")

N_per_dim=10
n_steps=10000

ρ=1.0
T=0.1

sys=sim_lennard_jones_fluid(N_per_dim,ρ,T,0.005,n_steps,VelocityVerlet,[],3.0,equilibration_steps=5000)
