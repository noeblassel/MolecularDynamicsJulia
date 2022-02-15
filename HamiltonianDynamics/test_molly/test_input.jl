using Molly
using Test,ProgressMeter

include("../utils/ReducedUnits.jl")
include("../molly/SimulateLennardJones.jl")
include("../molly/io.jl")
include("../molly/custom_simulators.jl")



ρ=0.5
T=0.5
N_per_dims=6
sys=sim_lennard_jones_fluid(N_per_dims,ρ,T,0.005,1000,VelocityVerlet,[],3.0;equilibration_steps=0)

save_reduced_lj_state(sys,"test_output.txt")
sys2=read_reduced_lj_state("test_output.txt")

@test length(sys)==length(sys2)

for i=1:length(sys)
    @test sys.coords[i]==sys2.coords[i]
    @test sys.velocities[i]==sys2.velocities[i]
end

simulate!(sys,SymplecticEulerA(dt=0.005),1000)
simulate!(sys2,SymplecticEulerA(dt=0.005),1000)

for i=1:length(sys)
    @test sys.coords[i]==sys2.coords[i]
    @test sys.velocities[i]==sys2.velocities[i]
end
