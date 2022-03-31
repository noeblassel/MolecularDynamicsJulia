using Plots

include("../molly/MollyExtend.jl")

Npd = 20
N = Npd^3
ρ = 1.0
L = (N / ρ)^(1 // 3)
T = 0.1
dt = 1e-3
r_c = 2.5

n_steps = 100_000_000

box_size = SVector(L, L, L)
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)

nf = CellListMapNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c,unit_cell=box_size)
coords =place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T, atoms[i].mass) for i in 1:N]


sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), box_size=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=Dict{Symbol,Any}())
sim_eq=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
sim=LangevinSplitting(dt=dt,γ=0.0,T=T,splitting="BABO")

simulate!(sys,sim_eq,20_000)

R(s::System,neighbors=nothing)=s.velocities[1][1]

sys.loggers=Dict(:autocorrelation=>AutoCorrelationLogger(R,1000))
simulate!(sys,sim,n_steps)
f=open("output.txt","w")
println(f,join(sys.loggers[:autocorrelation].correlations," "))
println(f,sys.loggers[:autocorrelation].avg_A," ",sys.loggers[:autocorrelation].avg_sq_A)
close(f)

"""
dt_range=1:dt:(dt*500)
scatter(dt_range,sys.loggers[:autocorrelation].correlations,color=:black,markershape=:circle,markersize=0.5,label="",xlabel="t",ylabel="C(t)",size=(1024,768))
savefig("solid_autocorrelation.png")
scatter(dt_range,cumsum(sys.loggers[:autocorrelation].correlations),color=:black,markershape=:circle,markersize=0.5,label="",xlabel="t",ylabel="∫C(t)dt",size=(1024,768))
savefig("int_solid_autocorrelation.png")"""