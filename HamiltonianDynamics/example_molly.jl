using Molly
using Unitful

include("simulation.jl")
include("animation.jl")

N_per_side = 6
T = 300u"K"#K
L = 7.0u"nm"#nm
ρ = (N_per_side^3 * m) / (L^3)
Δt = 0.001u"ps"
tfin=0.05u"ns"
obs = [("positions", 1)]
sys=sim_argon(N_per_side,ρ,T,Δt,tfin,VelocityVerlet,obs)


println(first(sys.loggers["positions"].coords))
scatter(first(sys.loggers["positions"].coords))
