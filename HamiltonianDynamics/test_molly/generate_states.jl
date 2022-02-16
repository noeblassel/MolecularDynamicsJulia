#!/libre/blasseln/julia-1.7.2/bin/julia

using Pkg

Pkg.add("Molly")
Pkg.add("ProgressMeter")

using Molly
using ProgressMeter

include("../molly/io.jl")
include("../molly/SimulateLennardJones.jl")

N_per_dim=10
ρ=1.0
Trange=0.01:0.1:10.0
n_steps=1000

for T in Trange
sys=sim_lennard_jones_fluid(N_per_dim,ρ,T,0.005,n_steps,VelocityVerlet,[],3.0,equilibration_steps=0)
save_reduced_lj_state(sys,"starting_states/config_T($(T)).lj_out")
end
