using Molly,ProgressMeter

include("../../molly/io.jl")
include("../../molly/sim_nve_lj.jl")
include("../../utils/reduced_units.jl")
include("../../utils/place_atoms.jl")

N_per_dim=10
ρ=1.0
Trange=0.01:0.1:10.0
n_steps=1000

for T in Trange
sys=sim_lennard_jones_fluid(N_per_dim,ρ,T,0.005,n_steps,VelocityVerlet,[],3.0,equilibration_steps=0)
save_reduced_lj_state(sys,"starting_states/config_T($(T)).lj_out")
end
