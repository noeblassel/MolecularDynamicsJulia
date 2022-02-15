
function save_reduced_lj_state(sim_result::NBodySimulator.SimulationResult,filename::AbstractString)

    @assert isa(sim_result.simulation.system,PotentialNBodySystem) "sim_result must contain the result of the simulation of a PotentialNBodySystem"
    @assert length(sim_result.simulation.system.potentials)==1 "The system can only contain more than one potential"
    @assert :lennard_jones in keys(sim_result.simulation.system.potentials) "The system contains a non Lennard-Jones potential"
    @assert isa(sim_result.simulation.boundary_conditions,CubicPeriodicBoundaryConditions) "Only supports cubic periodic boundary conditions"
    file=open(filename,"w")

    L=sim_result.simulation.boundary_conditions.L
    println(file,"$(L) $(L) $(L) $(sim_result.simulation.system.potentials[:lennard_jones].R)")

    t_fin=last(sim_result.simulation.tspan)

    N=length(sim_result.simulation.system.bodies)

    x=get_position(sim_result,t_fin)
    v=get_velocity(sim_result,t_fin)

    for i=0:N-1
        println(file,"$(x[3i+1]) $(x[3i+2]) $(x[3i+3]) $(v[3i+1]) $(v[3i+2]) $(v[3i+3])")
    end
    close(file)
end

function read_reduced_lj_state(filename::AbstractString)::NBodySimulation
    file=open(filename,"r")
    bodies=MassBody[]

    L=nothing
    R=nothing

    for (i,l) in enumerate(eachline(file))
        if i>1
            x1,x2,x3,v1,v2,v3=parse.(Float64,split(l))
            X=SVector(x1,x2,x3)
            V=SVector(v1,v2,v3)
            push!(bodies,MassBody(X,V,1.0))
        else
            L,_,_,R=parse.(Float64,split(l))
        end
    end

    potentials=Dict(:lennard_jones=>LennardJonesParameters(1.0,1.0,R))


    system=PotentialNBodySystem(bodies,potentials)
    boundary_conditions = CubicPeriodicBoundaryConditions(L)
    sim=NBodySimulation(system,(0.0,1.0),boundary_conditions,NBodySimulator.NullThermostat(),1.0)
    return sim
end

