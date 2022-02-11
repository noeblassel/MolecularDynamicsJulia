using Molly
using ProgressMeter

include("../molly/io.jl")
include("../molly/SimulateLennardJones.jl")

N_per_dim=6
n_steps=10000

ρ=0.5

for T in 0.01:0.01:0.5
    sys=sim_lennard_jones_fluid(N_per_dim,ρ,T,0.005,n_steps,VelocityVerlet,[],3.0,equilibration_steps=0)
    save_reduced_lj_state(sys,"starting_states/T($(round(T,digits=3))).out")
end

