using Molly
using Unitful
using GLMakie
using StaticArrays

Ar_mass=39.948u"u"
N=100
L=SVector(2.0,2.0,2.0)u"nm"
ρ=(N*Ar_mass)/(L[1]*L[2]*L[3])u"u * nm^-3"

particles=@SVector [Atom(mass=Ar_mass,σ=0.3u"nm",ϵ=0.2u"kJ * mol^-1") for i in 1:N]
coords=place_atoms(N,L,0.3u"nm")

T=120.0u"K"
velocities=@SVector [velocity(Ar_mass,T) for i in 1:N]

interactions=(LennardJones(),)
sys=System(atoms=particles,coords=coords,velocities=velocities,box_size=L,loggers=Dict("coords"=>CoordinateLogger(10),"total_energy"=>TotalEnergyLogger(10),"temperature"=>TemperatureLogger(10)))

simulator=StormerVerlet(dt=0.002u"ps")

simulate!(sys,simulator,1000)

visualize(sys.loggers["coords"],L,"traj.mp4")