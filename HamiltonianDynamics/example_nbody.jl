using NBodySimulator

include("nbodysimulator/SimulateLennardJones.jl")
include("utils/ReducedUnits.jl")

N_per_dim=10
T=0.5
ρ=0.5
Δt=0.005

n_steps=2000

res=sim_lennard_jones_fluid(N_per_dim,ρ,T,Δt,n_steps)
potential_energy