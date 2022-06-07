include("../molly/MollyExtend.jl")

using Plots, .MollyExtend

function ForcingProfileNorton(;n_bins::Integer,G::Function,dG::Function,γ::Float64,v::Float64,velocity_type::DataType=Float64)
    function R(sys::System,neighbors=nothing)
        bins=zeros(velocity_type,n_bins)
        N=length(sys)
        Ly=sys.box_size[2]
        accels=accelerations(sys,neighbors)
        for i=1:N
            Gy=G(sys.coords[i][2])
            bin_ix=1+floor(Int64,n_bins*ustrip(sys.coords[i][2]/Ly))
            bins[bin_ix]+=(v*dG(sys.coords[i][2])*sys.velocities[i][2]-accels[i][1]+γ*v*Gy)
        end
        return n_bins*bins/N
    end
    return R
end


G=ARGS[1]

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 16
N = Npd^3
γ=1.0
dt=5e-3
v=1.0
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

f_dict=Dict("SINUSOIDAL"=>(y-> sin(2π*y/L)),"CONSTANT"=>(y -> (y<L/2) ? -1 : 1),"LINEAR"=>(y -> (y<L/2) ? 4*(y-L/4)/L : 4*(3L/4-y)/L))
df_dict=Dict("SINUSOIDAL"=>(y-> (2π/L)*cos(2π*y/L)),"CONSTANT"=>(y -> 0),"LINEAR"=>(y -> (y<L/2) ? 4y/L : -4y/L))

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
simulator=NortonShearViscosityTest(dt = dt, γ = γ, T = T,v=v,G=f_dict[G])
loggers = Dict(:fp=>AverageObservableVecLogger(ForcingProfileNorton(n_bins=n_bins,G=f_dict[G],dG=df_dict[G],γ=γ,v=v),n_bins))

n_eq_steps=5000
n_steps=200000

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,),box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

simulate!(sys,simulator,n_eq_steps)

##reset loggers
sys.loggers[:fp].sum=zero(sys.loggers[:fp].sum)
sys.loggers[:fp].n_samples=0

println("equilibriated")

for i=1:20
    println(i)
    simulate!(sys,simulator,n_steps)
    fp_logger=sys.loggers[:fp]
    f_profile=fp_logger.sum/(fp_logger.n_samples)
    y_range=range(0,L,n_bins)
    plot(y_range,f_profile/v,label="",xlabel="y",ylabel="forcing",color=:red)
    plot!(f_dict[G],0,L,linestyle=:dot,color=:blue,label="")
    savefig("forcing_$(G)_norton.pdf")
end

