using Molly
using ProgressMeter
using Unitful, UnitfulRecipes
using Statistics
using Plots



include("molly/SimulateLennardJones.jl")
include("utils/ReducedUnits.jl")
include("utils/animate.jl")

N_per_side = 6
T = 0.71
ρ = 0.844
max_n = 1000
Δt = 0.005
obs = [(:state,500,"save_test")]
sys = sim_lennard_jones_fluid(N_per_side, ρ, T, Δt, max_n, SymplecticEulerA, obs, 3.0, equilibration_steps = 0)

