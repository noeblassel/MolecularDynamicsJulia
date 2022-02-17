using Molly,ProgressMeter,Plots,StaticArrays

include("../molly/custom_loggers.jl")
include("../molly/custom_cutoffs.jl")
include("../molly/custom_observables.jl")
configuration_file="../lj_sample_configurations/lj_sample_config_periodic1.txt"

r_c=3.0

file=open(configuration_file,"r")
(box,N,x...)=eachline(file)
coords=SVector{3}[]
velocities=SVector{3}[]
box_size=SVector(parse.(Float64,split(box))...)

for l ∈ x
    _,a,b,c=parse.(Float64,split(l))    
    push!(coords,SVector(a,b,c))
    push!(velocities,SVector(0.0,0.0,0.0))
end
N=parse(Int64,N)
println("$(N) atoms")
atoms=[Atom(mass=1.0,ϵ=1.0,σ=1.0) for i=1:N]
interactions = (LennardJones(cutoff = DistanceCutoff(r_c), force_units = NoUnits, energy_units = NoUnits,nl_only=false),)
nf=NoNeighborFinder()
nf=CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=1,dist_cutoff=r_c)
sys=System(atoms=atoms,general_inters=interactions,coords=coords,velocities=velocities,box_size=box_size,energy_units=NoUnits,force_units=NoUnits,neighbor_finder=nf)
neighbors=find_neighbors(sys,sys.neighbor_finder)

println(potential_energy(sys))  
println(pair_virial(sys,neighbors))
println(pressure(sys,neighbors))

println(box_size)