using Plots

include("../molly/MollyExtend.jl")

O(s::System,neighbors=nothing)=s.velocities[1][1]

Npd = 10
N = Npd^3
ρ = 0.4
L = (N / ρ)^(1 // 3)
T = 1.5
dt = 5e-3
r_c = 4.0

n_steps = 2000000

box_size = SVector(L, L, L)
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)

nf = CellListMapNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c,unit_cell=box_size)
coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T, atoms[i].mass) for i in 1:N]


sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), box_size=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=Dict{Symbol,Any}())
sim_eq=LangevinSplitting(dt=5e-3,γ=1.0,T=1.5,splitting="BAOAB")

simulate!(sys,sim_eq,5000)
sim=LangevinSplitting(dt=5e-3,γ=0.0,T=1.5,splitting="BAOAB")

sys.loggers=Dict(:autocorrelation=>TimeCorrelationLogger(O,O,5000))
simulate!(sys,sim,n_steps)

f=open("dump.out","r")
n=countlines(f)
seekstart(f)

dt_range=0:5e-3:(4999*5e-3)
anim=@animate for i=1:n
    C=parse.((Float64,),split(readline(f)))
    scatter(dt_range,C,markershape=:circle,color=:black,markersize=0.5,label="",ylims=(-0.2,1.2),size=(1024,768))
end

mp4(anim,"autocorrelations.mp4",fps=10)