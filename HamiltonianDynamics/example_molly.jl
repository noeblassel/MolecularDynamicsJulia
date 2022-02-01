using Molly
using Unitful,UnitfulRecipes,Plots
using Statistics
using GLMakie
"""
REFERENCE UNITS

mass: u
distance: nm
energy: kJ * mol ^-1

"""
##constants
const m = 39.948u"u"
const σ=0.3405u"nm"
const ϵ=1.0u"kJ * mol^-1"

function run_sim(N_per_dim,ρ,T,Δt,steps,log_coords=false)

    N = N_per_dim^3
    atoms = [Atom(mass=m, σ=σ, ϵ=ϵ) for i in 1:N]
    L=(N*m/ρ)^(1//3)
    l=L/N_per_dim
    domain = SVector(L, L, L)

    coords=[SVector(i*l,j*l,k*l) for i=1:N_per_dim,j=1:N_per_dim,k=1:N_per_dim]
    coords=reshape(coords,N)

    velocities = [velocity(m, T) for i in 1:N]
    interactions = (LennardJones(cutoff=DistanceCutoff(L/4)),)

    if log_coords
        loggers=Dict(
            "temp"   => TemperatureLogger(10),
            "total_energy" => TotalEnergyLogger(10),
            "coords"=>CoordinateLogger(100)
        )
    else
        loggers=Dict(
        "temp"   => TemperatureLogger(10),
        "total_energy" => TotalEnergyLogger(10),
    )
    end
    sys=System(atoms=atoms,general_inters=interactions,coords=coords,velocities=velocities,box_size=domain,loggers=loggers)

    simulator = VelocityVerlet(dt=Δt)
    simulate!(sys, simulator, steps),

    return sys
end

N_per_side=6
T = 100u"K"#K
L=50.0u"nm"#nm
ρ =(N_per_side^3*m)/(L^3)
Δt=0.0002u"ps"
sys=run_sim(N_per_side,ρ,T,Δt,1000000,true)

coords=sys.loggers["coords"].coords
X0=first(coords)
Xfin=last(coords)
println(X0-Xfin)

"""Δt_range=(0.0001u"ps":0.00001u"ps":0.002u"ps")

sim_duration=10000*0.0001u"ps"

energy_fluctuations=[]

for Δt=Δt_range
    n_steps=Int32(ceil(sim_duration/Δt))
    (H0,_,_,energies)=run_sim(N_per_side,ρ,T,Δt,n_steps)
    push!(energy_fluctuations,first(findmax(map(x->abs(x.val),H0.-energies))))
end

println(energy_fluctuations)
savefig(plot(Δt_range,energy_fluctuations),"energy_fluctuations.png")"""