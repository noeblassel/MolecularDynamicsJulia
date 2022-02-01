using Molly
using Unitful,UnitfulRecipes,Plots
using Statistics

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

function run_sim(N_per_dim,ρ,T,Δt,steps)

    N = N_per_dim^3
    atoms = [Atom(mass=m, σ=σ, ϵ=ϵ) for i in 1:N]
    L=(N*m/ρ)^(1//3)
    l=L/N_per_dim
    box_size = SVector(L, L, L)

    coords=[SVector(i*l,j*l,k*l) for i=1:N_per_dim,j=1:N_per_dim,k=1:N_per_dim]
    coords=reshape(coords,N)

    velocities = [velocity(m, T) for i in 1:N]
    interactions = (LennardJones(cutoff=DistanceCutoff(L/4)),)
    sys = System(
        atoms=atoms,
        general_inters=interactions,
        coords=coords,
        velocities=velocities,
        box_size=box_size,
        loggers=Dict(
            "temp"   => TemperatureLogger(1),
            "total_energy" => TotalEnergyLogger(1),
        ),
    )
    H0=total_energy(sys)
    time_range=(0.0u"ns":Δt:(steps-1)*Δt)
    simulator = VelocityVerlet(dt=Δt)
    simulate!(sys, simulator, steps)

    return (H0,time_range,sys.loggers["temp"].temperatures,sys.loggers["total_energy"].energies)
end

N_per_side=6
T = 100u"K"#K
L=50.0u"nm"#nm
ρ =(N_per_side^3*m)/(L^3)

Δt_range=(0.0001u"ps":0.00001u"ps":0.002u"ps")

sim_duration=10000*0.0001u"ps"

energy_fluctuations=[]

for Δt=Δt_range
    n_steps=Int32(ceil(sim_duration/Δt))
    (H0,_,_,energies)=run_sim(N_per_side,ρ,T,Δt,n_steps)
    push!(energy_fluctuations,first(findmax(map(x->abs(x.val),H0.-energies))))
end

println(energy_fluctuations)
savefig(plot(Δt_range,energy_fluctuations),"energy_fluctuations.png")