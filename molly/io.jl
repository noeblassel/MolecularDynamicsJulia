export 
    save_reduced_state,
    read_reduced_state

"""
Writes the state of a system to a given file-- only suitable for simple systems with general and pairwise interactions. For complex systems, use Molly's PDB writing utility.
"""
function save_reduced_state(s::System,file::IOStream)

    (length(s.general_inters)>0) && println(file,"---- General Interactions ----")
    for inter in values(s.general_inters)
        println(file,typeof(inter))
        dump(file,inter)
    end
    (length(s.pairwise_inters)>0) && println(file,"---- Pairwise Interactions -----")
    for inter in values(s.pairwise_inters)
        println(file,typeof(inter))
        dump(file,inter)
    end
    println(file,"---- Physical parameters ----")
    println(file,"N atoms : $(length(s))")
    println(file,"Box size : $(s.box_size)")
    
    println(file,"---- Atoms ----")
    for i=1:length(s)
        dump(file,s.atoms[i])
    end
    println(file,"---- Coordinates ----")
    for i = 1:length(s)
        x=s.coords[i]
        print(file,"[")
        println(file,"$(x[1]), $(x[2]), $(x[3]),")
    end
    print(file,"]")
    println(file,"---- Velocities ----")
    for i = 1:length(s)
        v=s.coords[i]
        print(file,"[")
        println(file,"$(v[1]), $(v[2]), $(v[3]),")
    end
    print(file,"]")
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

    if minimum(box_size)>3*r_c

        nf=CellListMapNeighborFinder(nb_matrix=trues(length(atoms),length(atoms)),dist_cutoff=r_c,unit_cell=box_size)
    else
        nf=TreeNeighborFinder(nb_matrix=trues(length(atoms),length(atoms)),dist_cutoff=r_c)
    end
    ##Assume ShiftedPotentialCutoff
    sys=System(atoms = atoms, pairwise_inters =(LennardJones(cutoff=ShiftedPotentialCutoff(r_c),force_units=NoUnits,energy_units=NoUnits),), coords = coords, velocities = velocities, box_size = box_size, energy_units = NoUnits, force_units = NoUnits)
    return sys
end

