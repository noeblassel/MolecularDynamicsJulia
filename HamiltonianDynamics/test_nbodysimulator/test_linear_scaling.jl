using Pkg

Pkg.add.(["NBodySimulator", "Unitful"])

using Unitful, NBodySimulator

include("../utils/ReducedUnits.jl")
include("../nbodysimulator/SimulateLennardJones.jl")

times = []

T = 1.0
ρ = 0.2
r_c = 3.0

Npd_range = 4:14
n_samps = 1
n_steps = 1000

for N = Npd_range
    println(N)
    t = 0
    for i = 1:n_samps
        t += @elapsed sim_lennard_jones_fluid(N, ρ, T, r_c, 5e-3, n_steps, VelocityVerlet)
    end
    push!(times, t / n_samps)
end

f = open("nbody_scaling.txt", "w")

for t in times
    print(f, "$(t) ")
end
println(f)

for Npd in 4:14
    print(f, "$(Npd^3) ")
end
close(f)
