using Molly, LinearAlgebra

include("../utils/place_atoms.jl")
println("Usage: T ρ dt γ η forcing_type=COLOR|SINGLE t_equilibration t_simulation N_atoms_per_dimension scheme cutoff_radius")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
η=parse(Float64,ARGS[5])
forcing_type=ARGS[6]
t_eq=parse(Float64,ARGS[7])
n_iter_sim=parse(Int64,ARGS[8])

Npd=parse(Int64,ARGS[9])
splitting=ARGS[10]
r_c=parse(Float64,ARGS[11])

N=Npd^3
L=(N/ρ)^(1//3)
box_size=CubicBoundary(L,L,L)

@assert 2r_c<= L "Cutoff radius too large relative to domain"


struct NEMD_longitudinal_forcing{F}
    forcing::F #transverse profile of the forcing
    η::Float64
end

function Molly.forces(inter::NEMD_longitudinal_forcing{F},s::System,neighbors=nothing) where {F}
    f=zero(s.velocities)
    f_x=view(reinterpret(reshape,Float64,f),1,:)
    q_y=view(reinterpret(reshape,Float64,s.coords),2,:)
    f_x .= inter.forcing.(q_y)
    return forcing.η*f
end

function fourier_response(s::System,args...;kwargs...)
    p_x=view(reinterpret(reshape,Float64,s.velocities),1,:)
    q_y=view(reinterpret(reshape,Float64,s.coords),2,:)
    Ly=s.boundary.side_lengths[2]
    N=length(s)
    return dot(p_x,exp.(2im*π*q_y/Ly))/N
end

sinus_forcing=(y-> sin(2π*y/L))
constant_forcing=(y -> (y<L/2) ? 1 : -1)
linear_forcing=(y -> (y<L/2) ? 4*(y-L/4)/L : 4*(3L/4-y)/L)

forcing=NEMD_longitudinal_forcing(sinus_forcing,η)

if forcing_type=="LINEAR"
    forcing=NEMD_longitudinal_forcing(linear_forcing,η)
elseif forcing_type=="CONSTANT"
    forcing=NEMD_longitudinal_forcing(constant_forcing,η)
end

nf = (3.6r_c < L) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),dist_cutoff=1.2r_c)

atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=[velocity(1.0,T,1.0) for i=1:N]
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))

sim=LangevinSplitting(dt=dt,friction=γ,temperature=T,splitting=splitting;remove_CM_motion=false)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)

simulate!(sys,sim,n_steps_eq)
sys=System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0,loggers=(fourier=GeneralObservableLogger(fourier_response,ComplexF64,1),))

for i=1:n_iter_sim
    simulate!(sys,sim,n_steps_eq)
    f=open("fourier_response_$(forcing_type)_$(η).out","a")
    write(f,values(sys.loggers.fourier))
    close(f)
    empty!(sys.loggers.fourier.history)
end