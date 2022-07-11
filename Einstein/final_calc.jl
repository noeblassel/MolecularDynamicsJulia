using Molly,LinearAlgebra

mutable struct SelfDiffusionLogger{T}
    log_every::Int64
    last_coords::T
    self_diffusion_coords::T
    msds::Vector{Float64}
end

SelfDiffusionLogger(initial_coords,log_every) = SelfDiffusionLogger{typeof(initial_coords)}(log_every,deepcopy(initial_coords), zero(initial_coords),Float64[])

function Molly.log_property!(logger::SelfDiffusionLogger, s::System, neighbors=nothing, step_n::Integer=0;parallel::Bool=true)
    logger.self_diffusion_coords += unwrap_coords_vec.(logger.last_coords, s.coords, (s.boundary.side_lengths,)) - logger.last_coords
    logger.last_coords .= s.coords
    (step_n % log_every == 0) && push!(logger.msds,ustrip(dot(logger.self_diffusion_coords,logger.self_diffusion_coords))/(3*length(logger.self_diffusion_coords)))
end

unwrap_coords(c1, c2, side_length)=c2+round((c1-c2)/side_length)*side_length
unwrap_coords_vec(c1, c2, box_size) = unwrap_coords.(c1, c2, box_size)

T=1.25
ρ=0.6
dt=1e-3
γ=1.0
t_eq=100.0
t_corr=2.0
Npd=10
splitting="BAOAB"
r_c=2.5
output_file="msds.out"

include("../utils/place_atoms.jl")

N=Npd^3
L=(N/ρ)^(1//3)
boundary=CubicBoundary(L,L,L)
@assert 2r_c<= L "Cutoff radius too large relative to domain"

nf = CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff= 1.2r_c,unit_cell=boundary)


atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,boundary)
velocities=[velocity(1.0,T,1.0) for i=1:N]
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))
l_corr=ceil(Int64,t_corr/dt)
n_iter_sim=20000

sim=LangevinSplitting(dt=dt,temperature=T,friction=γ,splitting="BAOAB",remove_CM_motion=false)
norm_cst=inv(sqrt(3N))
R(sys,args...;kwargs ...)=norm_cst*sys.velocities
sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),neighbor_finder=nf,boundary=boundary,energy_units=NoUnits,force_units=NoUnits,k=1.0)
simulate!(sys,sim,n_steps_eq)
sys=System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),neighbor_finder=nf,boundary=boundary,energy_units=NoUnits,force_units=NoUnits,k=1.0,loggers=(sd=SelfDiffusionLogger(sys.coords,100),))

for i=1:n_iter_sim
simulate!(sys,sim,n_steps_eq)
    f=open(output_file,"a")
    print(f,join(sys.loggers.sd.msds,'\n'))
    close(f)
    empty!(sys.loggers.sd.msds)
end

