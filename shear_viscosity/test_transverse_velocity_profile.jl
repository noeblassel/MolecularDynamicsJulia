include("../molly/MollyExtend.jl")

using Plots, .MollyExtend


F=ARGS[1]
ρ = 0.7
T = 1.0

r_c = 3.0

Ny = 15
ratio=5
N = Ny^3*ratio
γ=1.0
dt=5e-3
v=0.1
n_bins=300

function initialize_coords(ρ::Float64,Ny::Int64,ratio::Int64)
    L=Ny*(1/ρ)^(1//3)
    box_size=SVector(L*ratio,L,L)
    dy=L/Ny
    coords=dy*reshape([SVector{3,Float64}(x,y,z) for x=0:(Ny*ratio-1),y=0:(Ny-1),z=0:(Ny-1)],Ny^3*ratio)
    return coords,box_size
end

coords, box_size = initialize_coords(ρ,Ny,ratio)
Ly=box_size[2]
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = sqrt(T)*[SVector{3}(randn(3)) for i=1:N]

nf = nothing

if 3r_c<Ly
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

inter_dict=Dict("SINUSOIDAL"=>SinusoidalForceProfile,"LINEAR"=>PiecewiseLinearForceProfile,"CONSTANT"=>PiecewiseConstantForceProfile,"LOWER_BOUNDED"=>LowerBoundedSinusoidalForceProfile)

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
ff=inter_dict[F](ξ=ξ,L=L)
simulator=LangevinSplitting(dt = dt, γ = γ, T = T,splitting="BAOAB")
loggers = Dict(:vp=>AverageObservableVecLogger(TransverseVelocityProfile(n_bins=n_bins,L=L),n_bins))

n_eq_steps=5000
n_steps=100000

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), general_inters=(ff,),box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

simulate!(sys,simulator,n_eq_steps)

##reset loggers
sys.loggers[:vp].sum=zero(sys.loggers[:vp].sum)
sys.loggers[:vp].n_samples=0

println("equilibriated")

for i=1:100
    simulate!(sys,simulator,n_steps)
    vp_logger=sys.loggers[:vp]

    f=open("velocity_thevenin_$(F).out","w")
    println(f,"Ly: $Ly")
    println(f,"num_bins: $n_bins")
    println(f,"n_samples: $(vp_logger.n_samples)")
    println(f,"v: $v")
    println(f,join(vp_logger.sum," "))
    close(f)
end