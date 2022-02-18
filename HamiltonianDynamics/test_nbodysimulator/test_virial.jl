using NBodySimulator,Base.Threads


include("../nbodysimulator/custom_observables.jl")
include("../nbodysimulator/SimulateLennardJones.jl")

Npd=6
ρ=0.2
T=0.5

res=sim_lennard_jones_fluid(Npd,ρ,T,4.0,5e-3,100,VelocityVerlet)

@time V1=[pair_virial_lj(res,time) for time in res.solution.t]
@time V2=[pair_virial_lj(res,time,false) for time in res.solution.t]

@assert all(V1.≈V2)