function save_reduced_lj_state(s::System,filename::AbstractString)

    @assert length(s.general_inters)==1 "Trying to save a system with more than one GeneralInteraction"
    @assert isempty(s.specific_inter_lists) "Trying to save a system with specific interactions"
    @assert typeof(s.general_inters[1])<:LennardJones "Trying to save a non Lennard-Jones system"
    @assert typeof(s.general_inters[1].cutoff)<:ShiftedPotentialCutoff "Currently only implemented for Shifted potential cutoff"

    file=open(filename,"w")
    println(file,"$(s.box_size[1]) $(s.box_size[2]) $(s.box_size[3]) $(s.general_inters[1].cutoff.dist_cutoff)")
    for i = 1:length(s)
        x=s.coords[i]
        v=s.velocities[i]
        println(file,"$(x[1]) $(x[2]) $(x[3]) $(v[1]) $(v[2]) $(v[3])")
    end
    close(file)
end

function read_reduced_lj_state(filename::AbstractString)
    file=open(filename,"r")
    atoms=Atom{Float64,Float64,Float64,Float64}[]
    coords=SVector{3,Float64}[]
    velocities=SVector{3,Float64}[]
    box_size=nothing
    r_c=nothing

    for (i,l) in enumerate(eachline(file))
        if i>1
            x1,x2,x3,v1,v2,v3=parse.(Float64,split(l))
            push!(coords,SVector(x1,x2,x3))
            push!(velocities,SVector(v1,v2,v3))
            push!(atoms,Atom(index=i,mass=1.0,σ=1.0,ϵ=1.0))
        else
            l1,l2,l3,r=parse.(Float64,split(l))
            box_size=SVector{3,Float64}(l1,l2,l3)
            r_c=r
        end
    end
    close(file)
    ##Assume ShiftedPotentialCutoff
    sys=System(atoms = atoms, general_inters =(LennardJones(cutoff=ShiftedPotentialCutoff(r_c),force_units=NoUnits,energy_units=NoUnits),), coords = coords, velocities = velocities, box_size = box_size, energy_units = NoUnits, force_units = NoUnits)
    return sys
end

