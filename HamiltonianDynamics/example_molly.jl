using Molly
using Unitful

include("simulation.jl")
include("animation.jl")

N_per_side = 6
T = 300
L = 7.0
ρ = (N_per_side^3 * m) / (L^3)
Δt = 0.001
tfin=0.05
obs = [("positions", 1)]
sys=sim_argon(N_per_side,ρ,T,Δt,tfin,SymplecticEulerB,obs)
println(typeof(first(sys.loggers["positions"].coords)))