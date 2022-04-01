using Plots

include("../molly/MollyExtend.jl")

Npd = 7
N = Npd^3
ρ = 0.8
L = (N / ρ)^(1 // 3)
T = 0.1
dt = 5e-3
r_c = 2.5

n_steps = 1_000_000

box_size = SVector(L, L, L)
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)

nf = CellListMapNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c,unit_cell=box_size)

coords =place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
M=[a.mass for a=atoms]
velocities = init_velocities(T,M,1.0)

rseed=2022

sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), box_size=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=Dict{Symbol,Any}())
sim_eq=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
#sim=LangevinSplitting(dt=dt,γ=1.0,T=T,splitting="BAOAB")
sim=sim_eq
simulate!(sys,sim_eq,20_000)

R(s::System,neighbors=nothing)=s.coords/(L^3*length(s))

sys.loggers=Dict(:autocorrelation=>AutoCorrelationLoggerVec(N,3,R,10000))
simulate!(sys,sim,n_steps;log_progress=true,log_every=1000)
f=open("output.txt","w")
println(f,join(sys.loggers[:autocorrelation].correlations," "))
#println(f,sys.loggers[:autocorrelation].avg_A," ",sys.loggers[:autocorrelation].avg_sq_A)
close(f)


#plot(sys.loggers[:autocorrelation].correlations)
"""
dt_range=1:dt:(dt*500)
scatter(dt_range,sys.loggers[:autocorrelation].correlations,color=:black,markershape=:circle,markersize=0.5,label="",xlabel="t",ylabel="C(t)",size=(1024,768))
savefig("solid_autocorrelation.png")
scatter(dt_range,cumsum(sys.loggers[:autocorrelation].correlations),color=:black,markershape=:circle,markersize=0.5,label="",xlabel="t",ylabel="∫C(t)dt",size=(1024,768))
savefig("int_solid_autocorrelation.png")"""