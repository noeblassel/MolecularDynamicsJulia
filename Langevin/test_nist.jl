using Statistics

include("../molly/MollyExtend.jl")

if length(ARGS)!=9
    println("Error parsing arguments.")
    println("Usage: julia test_nist.jl TEMPERATURE ρMIN ρMAX NUMBER_OF_SIMULATIONS ATOMS_PER_DIM EQUILIBRATION_STEPS SAMPLING_STEPS Δt FILEOUT")
    exit(1)
end
    
T = parse(Float64,ARGS[1])
ρmin=parse(Float64,ARGS[2])#0.4574955528339076
ρmax=parse(Float64,ARGS[3])#0.6858732065373103
n_sims=parse(Int32,ARGS[4])#30
Npd=parse(Int32, ARGS[5])
eq_steps=parse(Int32,ARGS[6])
samp_steps=parse(Int32, ARGS[7])
dt = parse(Float64,ARGS[8])
filename=ARGS[9]


ρs=ρmin:((ρmax-ρmin)/n_sims):ρmax


N = Npd^3
r_c=5.0

atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

for ρ in ρs
    
    L = (N / ρ)^(1 // 3)
    
    box_size = SVector(L, L, L)
    
    coords = place_atoms_on_lattice(Npd, box_size)
    
    velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]
    
    inter = LennardJones(cutoff = DistanceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c,unit_cell=box_size)   
    
    γ = 1.0
 
    
    simulator = LangevinBAOAB(dt = dt,γ=γ,T=T)
    
    sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)
    try 
        simulate!(sys, simulator, eq_steps)
        
        loggers = Dict(:pressure => PressureLoggerNVT(T,Float64, 1))

        sys = System(atoms = sys.atoms, coords = sys.coords, velocities = sys.velocities, pairwise_inters = (inter,), box_size = sys.box_size, neighbor_finder = sys.neighbor_finder, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)
        
        simulate!(sys, simulator, samp_steps)

        if filename!="STDOUT"
        f=open(filename,"a")
        println(f,ρ," ", mean(sys.loggers[:pressure].pressures)," ",long_range_virial_correction(sys,sys.pairwise_inters[1])/(3L^3))
        close(f)
        else
            println(ρ," ", mean(sys.loggers[:pressure].pressures)," ",long_range_virial_correction(sys,sys.pairwise_inters[1])/(3L^3))
        end
    catch
        if filename!="STDOUT"
        f=open(filename,"a")
        println(f,ρ," :error")
        close(f)
        else
            println(ρ," :error")
        end
    end
end
