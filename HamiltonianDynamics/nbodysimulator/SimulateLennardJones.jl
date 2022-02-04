using NBodySimulator
using StaticArrays

include("../ReducedUnits.jl")

function simulate_lennard_jones_fluid(species,N_per_dim, ρ, T, Δt, steps)
    N = N_per_dim^3
    L = (N * m * ρ)^(1.0 / 3)
    
    domain = SVector(L, L, L)
    v=sqrt(kb*T/m)
    bodies = generate_bodies_in_cell_nodes(N,m, v,L/N_per_dim)
    lj_params=LennardJonesParameters(1.0,1.0,2.5)

    system = PotentialNBodySystem(bodies,Dict(:lennard_jones=>lj_params))
    boundary_conditions=PeriodicBoundaryConditions(L)
    sim=NBodySimulation(system,(0.0,1.0),boundary_conditions)
    result=run_simulation(sim, VelocityVerlet(),dt=Δt)
    return result

end



#test run
N=5^3
L=2.0#nm
T = 200#K
ρ =L^3/(N*m)
Δt = 0.02#ns
res = run_sim(5,ρ,T,Δt,1000,L/5)
print(res)