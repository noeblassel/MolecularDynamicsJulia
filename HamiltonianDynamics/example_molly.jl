using Molly
using ProgressMeter
using Unitful,UnitfulRecipes
using Statistics
using Plots



include("molly/SimulateLennardJones.jl")
include("utils/ReducedUnits.jl")
include("utils/animate.jl")   

N_per_side = 8
T=0.71
ρ=0.844
max_n=10000
Δt=0.005
obs=[(:temperature,1),(:potential_energy,1),(:pressure,1)]
sys=sim_lennard_jones_fluid(N_per_side,ρ,T,Δt,max_n,VelocityVerlet,obs,3.0,equilibration_steps=10000)

(bins,dens)=Molly.rdf(sys.coords,sys.box_size)
plot(bins,dens)
