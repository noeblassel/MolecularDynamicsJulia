include("../molly/MollyExtend.jl")

using Plots, .MollyExtend

Npd=16

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3
γ=1.0
dt=5e-3
ξ=1.0
n_bins=100

L = (N / ρ)^(1 // 3)
box_size = SVector(L, L, L)

coords = place_atoms_on_3D_lattice(Npd,box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = init_velocities(T,[a.mass for a=atoms],1.0)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
ff=SinusoidalForceProfile(ξ=ξ,L=L)
simulator=LangevinSplitting(dt = dt, γ = γ, T = T,splitting="BAOAB")
loggers = Dict(:vp=>AverageObservableVecLogger(TransverseVelocityProfile(n_bins=n_bins,L=L),n_bins))

n_eq_steps=5000
n_steps=20000

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), general_inters=(ff,),box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

simulate!(sys,simulator,n_eq_steps)

##reset loggers
sys.loggers[:vp].sum=zero(sys.loggers[:vp].sum)
sys.loggers[:vp].n_samples=0

println("equilibriated")

for i=1:20
    println(i)
    simulate!(sys,simulator,n_eq_steps)
    vp_logger=sys.loggers[:vp]
    profile=vp_logger.sum/(vp_logger.n_samples*ξ)
    y_range=range(0,L,n_bins)
    plot(profile,y_range,label="computed profile",xlabel="x velocity",ylabel="y coordinate",color=:red)
    plot!(sin.(2π*y_range/L),y_range,label="force profile",linestyle=:dot,color=:blue)
    savefig("sinusoidal$i.pdf")
end