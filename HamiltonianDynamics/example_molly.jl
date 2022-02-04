using Molly
using Unitful

include("molly/SimulateLennardJones.jl")

N_per_side = 6
T = 60u"K"
L = 6.0u"nm"
ρ = (N_per_side^3)/(L^3)
Δt = 0.001
n_steps=100000
obs = [(:potential_energy,1),(:hamiltonian,1),(:kinetic_energy,1)]

sys=sim_lennard_jones_fluid(:Ar,N_per_side,ρ,T,Δt,n_steps,SymplecticEulerA,obs)

using Plots

plot(sys.loggers[:potential_energy].energies)
plot!(sys.loggers[:kinetic_energy].energies)
plot!(sys.loggers[:hamiltonian].energies)
savefig("equilibration.png")
