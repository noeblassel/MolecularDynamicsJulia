using Molly
using Unitful, UnitfulRecipes, Plots
using Statistics

include("simulation.jl")

N_per_side = 6
T = 300u"K"#K
L = 50.0u"nm"#nm
ρ = (N_per_side^3 * m) / (L^3)
Δt = 0.0002u"ps"
obs = [("hamiltonian", 10)]

Δt_range = (0.0001u"ps":0.0001u"ps":0.002u"ps")
println(length(Δt_range))

sim_duration = 0.1u"ns"
energy_fluctuations = []

for Δt = Δt_range
    n_steps = Int32(ceil(sim_duration / Δt))
    s = sim_argon(N_per_side, ρ, T, Δt, sim_duration,VelocityVerlet, obs)
    energy_record = s.loggers["hamiltonian"].energies
    time_range=(0.0u"ps":Δt:Δt*(n_steps-1))
    savefig(plot(energy_record), "./plots/$(Δt.val).png")
end
