using NBodySimulator, BenchmarkTools, Plots

include("../utils/ReducedUnits.jl")
include("../nbodysimulator/SimulateLennardJones.jl")
include("../nbodysimulator/io.jl")

times = []

T = 1.0
ρ = 0.2
r_c = 3

steps = 5000

Npd_range = 4:15
n_samps = 10

for N = 1:12
    t = 0
    for i = 1:n_samps
        t += @elapsed sim_lennard_jones_fluid(N, ρ, T, r_c, 5e-3, 10000, VelocityVerlet)
    end
    push!(times, t / n_samps)
end

f = open("nbody_scaling.txt", "w")

for t in times
    print(f, "$(t) ")
end
close(f)