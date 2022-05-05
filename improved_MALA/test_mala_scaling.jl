using Molly,LinearAlgebra,Random

include("../molly/MollyExtend.jl")
include("integrators.jl")

using .MollyExtend

println(stderr,"Usage: Npd log_dt_min log_dt_max num_dts ρ N_steps Nsamps rule(METROPOLIS|BARKER) proposal(EM|HMC)")
Npd=parse(Int64,ARGS[1])
lg_dt_min=parse(Float64,ARGS[2])
lg_dt_max=parse(Float64,ARGS[3])
N_dts=parse(Int64,ARGS[4])
ρ=parse(Float64,ARGS[5])
N_steps=parse(Int64,ARGS[6])
Nsamps=parse(Int64,ARGS[7])
rule=ARGS[8]
proposal=ARGS[9]

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
dt_eq=1e-3
sim_eq=MALA(dt=dt_eq, T=1.0)
simulate!(sys,sim_eq,n_steps_eq)
#---------- simulate --------------

metropolis= (rule=="METROPOLIS")

log_dts=range(lg_dt_min,lg_dt_max,N_dts)
dts= 10 .^ log_dts
As=zero(dts)
Rs_baker_abs=zero(dts)

for it=1:Nsamps
    println("iteration $it/$Nsamps")
    flush(stdout)
    for (i,dt)=enumerate(dts)

        if proposal=="EM"
            sim=MALA(dt=dt,T=1.0,is_metropolis=metropolis)
        else
            sim=MALA_HMC(dt=dt,T=1.0,is_metropolis=metropolis)
        end

        simulate!(sys,sim,N_steps)
        global As[i]+=(sim.n_accepted/sim.n_total)
        global Rs_baker_abs[i]+=sim.baker_abs_sum
    end
end

println(dts)

if metropolis
    println(1 .- As/Nsamps)

else
    println(abs.(2*As/Nsamps .- 1))
    println(Rs_baker_abs/(Nsamps*N_steps))
end