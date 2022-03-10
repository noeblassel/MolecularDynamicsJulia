include("../molly/MollyExtend.jl")
using Statistics

#julia test_bias.jl T ρ dt tfin simulator output

@assert length(ARGS)==6 "Bad argcount. Usage: test_bias.jl T ρ Δt Tfin BAOAB|BABO|BAOA OUTPUT|STDOUT"
T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
tfin=parse(Float64,ARGS[4])
sim=ARGS[5]
output_file=ARGS[6]

sys=read_reduced_lj_state("bias_preconfig/initial_config_T_$(T)_rho_$(ρ).out")


n_steps=Int64(round(tfin/dt))
println(n_steps)

loggers=Dict(:potential_energy=>PotentialEnergyLogger(Float64,1),:kinetic_energy=>KineticEnergyLoggerNoDims(Float64,1),:virial=>VirialLogger(Float64,1))

sys.loggers=loggers
pairwise_inters=(LennardJones(cutoff = ShiftedForceCutoff(sys.pairwise_inters[1].cutoff.dist_cutoff), nl_only = true, force_units = NoUnits, energy_units = NoUnits),)
nf=DistanceNeighborFinder(dist_cutoff=pairwise_inters[1].cutoff.dist_cutoff)


sys=System(atoms=sys.atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=pairwise_inters,loggers=sys.loggers,neighbor_finder=nf,box_size=sys.box_size,energy_units=NoUnits,force_units=NoUnits)
γ=1.0
simulator=LangevinBAOAB(dt=dt,T=T,γ=γ,rseed=123)


if sim=="BAOAB"
    simulator=LangevinBAOAB(dt=dt,T=T,γ=γ)
elseif sim=="BABO"
    simulator=LangevinBABO(dt=dt,T=T,γ=γ)
elseif sim=="BAOA"
    k=0.008314462621026539#to deal with Molly's assumption about units
    simulator=LangevinBAOA(dt=dt,temperature=T/k,friction=γ)
else
    println("unrecognized simulator.")
    exit(1)
end

simulate!(sys,simulator,n_steps)

f=open(output_file,"a")
Vhat=mean(sys.loggers[:potential_energy].energies)
Khat=mean(sys.loggers[:kinetic_energy].energies)
What=mean(sys.loggers[:virial].energies)

println(f,"$(dt) $(Vhat) $(Khat) $(What)")
close(f)