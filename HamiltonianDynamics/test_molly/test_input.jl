using Molly
using Test,ProgressMeter

include("../../utils/reduced_units.jl")
include("../../molly/sim_nve_lj.jl")
include("../../molly/io.jl")
include("../../molly/custom_simulators.jl")



ρ=0.5
T=0.5
N_per_dims=10
sys_a=sim_lennard_jones_fluid(N_per_dims,ρ,T,0.005,1000,VelocityVerlet,[],3.0;equilibration_steps=0)

save_reduced_lj_state(sys_a,"test_output.txt")
sys_b=read_reduced_lj_state("test_output.txt")

@test length(sys_a)==length(sys_b)

for i=1:length(sys_a)
    @test sys_a.coords[i]==sys_b.coords[i]
    @test sys_a.velocities[i]==sys_b.velocities[i]
end

@time simulate!(sys_a,VelocityVerlet(dt=0.005),1000)
@time simulate!(sys_b,VelocityVerlet(dt=0.005),1000)

for i=1:length(sys_a)
    @test sys_a.coords[i]==sys_b.coords[i]
    @test sys_a.velocities[i]==sys_b.velocities[i]
end
