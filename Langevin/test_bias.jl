include("../molly/MollyExtend.jl")
using Statistics

#julia test_bias.jl T ρ dt tfin simulator output

@assert length(ARGS)==10 "Error (Wrong Argument Count) Usage: test_bias.jl T ρ Δt teq tfin Nruns Npd r_c BAOAB|BABO|BAOA OUTPUT|STDOUT"
T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
teq=parse(Float64,ARGS[4])
tfin=parse(Float64,ARGS[5])
Npd=parse(Int64,ARGS[6])
Nruns=parse(Int64,ARGS[7])
r_c=parse(Float64,ARGS[8])
sim=ARGS[9]
output_file=ARGS[10]

if !isfile(output_file)
    f=open(output_file,"w")
    println(f,"Trajectorial averages of NVT Lennard-Jones system of $(Npd^3) particles at T=$(T), ρ=$(ρ) with $(sim) splitting, $(r_c) shifted force cutoff. Physical time of each run: $(tfin). All units are reduced.")
    println(f,"[dt] [average potential energy] [average kinetic energy] [average virial]")
    close(f)
end

dt_eq=5e-3
eq_nsteps=Int64(round(teq/dt_eq))
N = Npd^3
L = (N / ρ)^(1 // 3)
box_size = SVector(L, L, L)
inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf=nothing
n_steps=Int64(round(tfin/dt))

if L>3*r_c
    nf=CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c,unit_cell=box_size)
elseif N>900
    nf=TreeNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c)
else
    nf=DistanceNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c)
end

for i=1:Nruns
    #equilibriate
    coords = place_atoms_on_lattice(Npd, box_size)
    atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
    velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]


    sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits,loggers=Dict{Symbol,Any}())


    γ=1.0
    simulator=LangevinBAOAB(T=T,γ=γ,dt=dt_eq)
    simulate!(sys,simulator,eq_nsteps)

    loggers=Dict(:potential_energy=>PotentialEnergyLogger(Float64,1),:kinetic_energy=>KineticEnergyLoggerNoDims(Float64,1),:virial=>VirialLogger(Float64,1))

    sys.loggers=loggers
    simulator=LangevinBAOAB(dt=dt,T=T,γ=γ)

    if sim=="BAOAB"
        simulator=LangevinBAOAB(dt=dt,T=T,γ=γ)
    elseif sim=="BABO"
        simulator=LangevinBABO(dt=dt,T=T,γ=γ)
    elseif sim=="BAOA"
        k=0.008314462621026539#to deal with Molly's assumption about units
        simulator=LangevinBAOA(dt=dt,temperature=T/k,friction=γ)
    else
        println("unrecognized simulator.")
        exit(1)
    end

    simulate!(sys,simulator,n_steps)

    f=open(output_file,"a")
    Vhat=mean(sys.loggers[:potential_energy].energies)
    Khat=mean(sys.loggers[:kinetic_energy].energies)
    What=mean(sys.loggers[:virial].energies)

    println(f,"$(dt) $(Vhat) $(Khat) $(What)")
    close(f)
end