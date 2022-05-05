using Molly,LinearAlgebra,Random

include("../molly/MollyExtend.jl")
include("integrators.jl")

using .MollyExtend

println(stderr,"Usage: Npd dt ρ N_steps T_corr rule(METROPOLIS|BARKER) proposal(EM|HMC)")
Npd=parse(Int64,ARGS[1])
lg_dt=parse(Float64,ARGS[2])
ρ=parse(Float64,ARGS[3])
N_steps=parse(Int64,ARGS[4])
T_corr=parse(Float64,ARGS[5])
rule=ARGS[6]
proposal=ARGS[7]

N=Npd^3
dt=10^lg_dt
N_corr=ceil(Int64,T_corr/dt)
(rule == "BARKER") && (N_corr*=2)
R(sys::System,neighbors=nothing)=forces(sys,neighbors)/N #observable for autocorrelation,normalized by N for numerical stability

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

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c_inter), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
atoms=[Atom(ϵ=1.0,mass=1.0,σ=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Npd,box_size)
velocities=zero(coords)#not needed for overdamped langevin/ abused in the HMC proposal case
sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = Dict{Symbol,Any}())
#---------------------------equilibriate---------------------------------------

n_steps_eq=10_000
dt_eq=5e-4
sim_eq=MALA(dt=dt_eq, T=1.0)
simulate!(sys,sim_eq,n_steps_eq)
#---------- simulate --------------
sys.loggers=Dict(:autocorr=>AutoCorrelationLoggerVec(N,3,R,N_corr),:msd=>SelfDiffusionLogger(coords),:log=>LogLogger([:autocorr,:msd],["output/autocorrelation_$(rule)_$(proposal)_$(lg_dt).out","output/sd_increments_$(rule)_$(proposal)_$(lg_dt)"],[1000*N_corr,1000*N_corr],[false,false],["a","a"]))
metropolis = (rule=="METROPOLIS")
N_sim_steps=N_steps*N_corr
for i in 1 #to have access to sim in local scope
    if proposal=="EM"
        sim=MALA(dt=dt,T=1.0,is_metropolis=metropolis)
    else
        sim=MALA_HMC(dt=dt,T=1.0,is_metropolis=metropolis)
    end
    simulate!(sys,sim,N_sim_steps)
end