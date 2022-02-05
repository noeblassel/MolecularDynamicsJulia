using Molly
using Unitful
using LaTeXStrings

include("molly/SimulateLennardJones.jl")
include("ReducedUnits.jl")

N_per_side = 10
T=0.5
ρ=0.5
Δt = 0.02
n_steps=2500
obs = [(:position,1)]

sys=sim_lennard_jones_fluid(N_per_side,ρ,T,Δt,n_steps,SymplecticEulerA,obs)

include("animate.jl")

animate_trajectories(sys,"lj_test_symplecticEulerA_big_dt.gif")
#savefig(plot(sys.loggers[:temperature].temperatures),"kinetic_temperature.png")
