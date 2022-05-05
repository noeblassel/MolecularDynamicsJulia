using Molly,LinearAlgebra,Random

include("../molly/MollyExtend.jl")
include("integrators.jl")

using .MollyExtend

println(stderr,"Usage: Npd log_dt_min log_dt_max num_dts ρ N_steps Nsamps is_metropolis")
Npd=parse(Int64,ARGS[1])
lg_dt_min=parse(Float64,ARGS[2])
lg_dt_max=parse(Float64,ARGS[3])
N_dts=parse(Int64,ARGS[4])
ρ=parse(Float64,ARGS[5])
N_steps=parse(Int64,ARGS[6])
Nsamps=parse(Int64,ARGS[7])
metropolis=parse(Bool,ARGS[8])

N=Npd^3

#------------------------ setup neighbor list ----------------------------------
L = (N / ρ)^(1 // 3)
box_size=SVector(L,L,L)
r_c_nf=1.8

if (L > 3 * r_c_nf) && N>900
    nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf, unit_cell = box_size)
elseif N > 900
    nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf)
else
    nf = DistanceNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c_nf)
end
#----------------------- setup system -----------------------------------
r_c_inter=1.6

inter = SoftSphere(cutoff = ShiftedForceCutoff(r_c_inter), nl_only = true, force_units = NoUnits, energy_units = NoUnits)


atoms=[Atom(ϵ=1.0,mass=1.0,σ=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=init_velocities(1.0,[a.mass for a=atoms],1.0)
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = Dict{Symbol,Any}())
#---------------------------equilibriate---------------------------------------

n_steps_eq=10_000
dt_eq=5e-3
sim_eq=MALA_HMC(dt=dt_eq, T=1.0)
simulate!(sys,sim_eq,n_steps_eq)
#---------- simulate --------------

log_dts=range(lg_dt_min,lg_dt_max,N_dts)
dts= 10 .^ log_dts
Rs=zero(dts)
for it=1:Nsamps
    println("iteration $it/$Nsamps")
    for (i,dt)=enumerate(dts)
        sim=MALA_HMC(dt=dt,T=1.0,is_metropolis=metropolis)
        simulate!(sys,sim,N_steps)
        global Rs[i]+=1-(sim.n_accepted/sim.n_total)
    end
end

println(dts)
println(Rs/Nsamps)

julia test_mala_scaling.jl 4 -6 -3 20 0.4 100000 10 true