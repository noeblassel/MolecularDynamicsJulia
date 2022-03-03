include("../../molly/MollyExtend.jl")

N_per_dim=10
ρ=0.4
Trange=[1.0]
n_steps=10000

for T in Trange
sys=sim_lennard_jones_fluid_nve(N_per_dim,ρ,T,0.005,n_steps,VelocityVerlet,[],4.0,equilibration_steps=0)
save_reduced_lj_state(sys,"starting_states/config.lj_out")
end
