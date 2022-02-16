using Molly,ProgressMeter,Plots,StaticArrays

include("../molly/custom_loggers.jl")

configuration_file="../lj_sample_configurations/lj_sample_config_periodic1.txt"


file=open(configuration_file,"r")
(box,N,x...)=eachline(file)
coords=[]
velocities=[]
box_size=SVector(parse.(Float64,split(box))...)

for l ∈ x
    _,a,b,c=parse.(Float64,split(l))    
    push!(coords,SVector(a,b,c))
    push!(velocities,SVector(0.0,0.0,0.0))
end
N=parse(Int64,N)
atoms=[Atom(mass=1.0,ϵ=1.0,σ=1.0) for i=1:N]
interactions = (LennardJones(cutoff = ShiftedForceCutoff(3.0), force_units = NoUnits, energy_units = NoUnits),)
sys=System(atoms=atoms,general_inters=interactions,coords=coords,velocities=velocities,box_size=box_size,energy_units=NoUnits,force_units=NoUnits)

println(potential_energy(sys))  
println(virial(sys))