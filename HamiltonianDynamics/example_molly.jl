using Molly
using Unitful,UnitfulRecipes
using Statistics
using Plots



include("molly/SimulateLennardJones.jl")
include("ReducedUnits.jl")
include("animate.jl")

N_per_side = 10
T=0.1
ρ=0.9
Δt = 0.005
n_steps=2000
obs = [(:position,1),(:pressure,50),(:temperature,50),(:kinetic_energy,50),(:potential_energy,50)]
sys=sim_lennard_jones_fluid(N_per_side,ρ,T,Δt,n_steps,VelocityVerlet,obs)

(bins,densities)=rdf(sys.coords,sys.box_size)
savefig(plot(bins,densities),"test_rdf.png")


