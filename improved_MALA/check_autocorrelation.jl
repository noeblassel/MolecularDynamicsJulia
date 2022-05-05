using Molly,LinearAlgebra,Random

include("../molly/MollyExtend.jl")
include("integrators.jl")

using .MollyExtend

println(stderr,"Usage: Npd dt ρ N_steps T_corr rule(METROPOLIS|BAKER) proposal(EM|HMC)")
Npd=parse(Int64,ARGS[1])
dt=parse(Float64,ARGS[2])
ρ=parse(Float64,ARGS[3])
N_steps=parse(Int64,ARGS[4])
T_corr=parse(Float64,ARGS[5])
rule=ARGS[6]
proposal=ARGS[7]

N=Npd^3

N_steps_corr=ceil(Int64,T_corr/dt)
R(sys::System,neighbors=nothing)=forces(sys,neighbors)#observable for autocorrelation

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
sys.loggers=Dict(:autocorr=>AutoCorrelationLoggerVec(N,3,R,N_steps_corr),:log=>LogLogger([:autocorr],["autocorrelation_$rule_$proposal.out"],[1000000],[false],["w"]))
metropolis= (rule=="METROPOLIS")

for i in 1
    if proposal=="EM"
        sim=MALA(dt=dt,T=1.0,is_metropolis=metropolis)
    else
        sim=MALA_HMC(dt=dt,T=1.0,is_metropolis=metropolis)
    end
    simulate!(sys,sim,N_steps)
end