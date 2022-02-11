function read_reduced_lj_state(filename::AbstractString)
    file=open(filename,"r")
    atoms=Atom{Float64,Float64,Float64,Float64}[]
    coords=SVector{3,Float64}[]
    velocities=SVector{3,Float64}[]
    box_size=nothing
    for (i,l) in enumerate(eachline(file))
        if i>1
            x1,x2,x3,v1,v2,v3=parse.(Float64,split(l))
            push!(coords,SVector(x1,x2,x3))
            push!(velocities,SVector(v1,v2,v3))
            push!(atoms,Atom(index=i,mass=1.0,σ=1.0,ϵ=1.0))
        else
            box_size=SVector{3,Float64}(parse.(Float64,split(l)))
        end
    end
    close(file)
    ##Assume ShiftedPotentialCutoff, r_c=3.0
    sys=System(atoms = atoms, general_inters =(LennardJones(cutoff=ShiftedPotentialCutoff(3.0),force_units=NoUnits,energy_units=NoUnits),), coords = coords, velocities = velocities, box_size = box_size, energy_units = NoUnits, force_units = NoUnits)
    return sys
end
