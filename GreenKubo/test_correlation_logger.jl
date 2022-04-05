include("../molly/MollyExtend.jl")
using .MollyExtend

Npd = 12
N = Npd^3
ρ = 0.8;
L = (N / ρ)^(1 // 3)
T = 0.2
dt = 2e-3
r_c = 3.0

n_steps = 5_000_000

box_size = SVector(L, L, L)
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)
nf = CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=r_c,unit_cell=box_size)
coords =place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
M=[a.mass for a=atoms]
velocities = init_velocities(T,M,1.0)

sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), box_size=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=Dict{Symbol,Any}())
sim_eq=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
#sim=LangevinSplitting(dt=dt,γ=0.1,T=T,splitting="BAOAB")
sim=sim_eq
simulate!(sys,sim_eq,50_000)

R(s::System,neighbors=nothing)=s.velocities/N

sys.loggers=Dict(:autocorrelation=>AutoCorrelationLoggerVec(N,3,R,2000),:time=>ElapsedTimeLogger(),:state=>StateLogger(sys),:log_log=>LogLogger([:autocorrelation,:time,:state],["autocorrelation_history.out","elapsed_times.out","states.out"],[100_000,1000,1000000],[false,true]))

simulate!(sys,sim,n_steps)
f=open("output_alt.txt","w")
println(f,join(sys.loggers[:autocorrelation].correlations," "))
close(f)

##clustern19 :dt=2e-3
##clustern24 :dt=1e-3