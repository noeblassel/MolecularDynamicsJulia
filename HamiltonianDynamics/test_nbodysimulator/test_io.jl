using NBodySimulator,StaticArrays,Test

include("../nbodysimulator/io.jl")
include("../nbodysimulator/SimulateLennardJones.jl")
include("../utils/ReducedUnits.jl")

N_pd=6

ρ=0.5
T=0.5

res=sim_lennard_jones_fluid(N_pd,ρ,T,0.005,200,VelocityVerlet)

save_reduced_lj_state(res,"ini_config.out")
sim=read_reduced_lj_state("ini_config.out")

res1=sim_lennard_jones_fluid(sim,0.005,200,VelocityVerlet)
save_reduced_lj_state(res1,"result1.out")

res_a=sim_lennard_jones_fluid(sim,0.005,400,VelocityVerlet)
sim2=read_reduced_lj_state("result1.out")

res_b=sim_lennard_jones_fluid(sim2,0.005,200,VelocityVerlet)

@test NBodySimulator.get_coordinate_vector_count(res_a.simulation.system)==NBodySimulator.get_coordinate_vector_count(res_a.simulation.system)
N=NBodySimulator.get_coordinate_vector_count(res_a.simulation.system)
t_f_a=last(res_a.simulation.tspan)
t_f_b=last(res_b.simulation.tspan)

x_a=get_position(res_a,t_f_a)
x_b=get_position(res_b,t_f_b)

v_a=get_velocity(res_a,t_f_a)
v_b=get_velocity(res_b,t_f_b)

@test x_a ≈ x_b atol=1e-5
@test v_a ≈ v_b atol=1e-5

