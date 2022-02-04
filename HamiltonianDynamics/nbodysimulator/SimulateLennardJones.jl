using NBodySimulator
using StaticArrays


#
const kb = 8.3144598e-3
const Nₐ = 6.0221415*1023

#Lennard Jones parameters for Argon
const σ = 0.3405#nm
const ϵ = 1.66e-21 * Nₐ#J* mol^-1

#mass of an Argon atom
const m = 6.634e-26 #kg

function run_sim(N_per_dim, ρ, T, Δt, steps, r_cutoff)
    N = N_per_dim^3
    L = (N * m * ρ)^(1.0 / 3)
    
    domain = SVector(L, L, L)
    v=sqrt(kb*T/m)
    bodies = generate_bodies_in_cell_nodes(N,m, v,L/N_per_dim)
    lj_params=LennardJonesParameters(ϵ,σ,r_cutoff)

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