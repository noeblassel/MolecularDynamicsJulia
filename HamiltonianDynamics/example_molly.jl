using Molly
using Unitful
using LaTeXStrings

include("molly/SimulateLennardJones.jl")
include("ReducedUnits.jl")

N_per_side = 6
T=1.5043
ρ=0.8442
Δt = 0.003
n_steps=2500
obs = [(:position,1)]

sys=sim_lennard_jones_fluid(N_per_side,ρ,T,Δt,n_steps,VelocityVerlet,obs)

include("animate.jl")

animate(sys,"lj_test.gif")
