include("../molly/MollyExtend.jl")
using .MollyExtend

Npd=3
rseed=2022
N=Npd^3
T=1.0

bs_1=SVector(2.0,2.0,2.0)
bs_2=SVector(1.0,1.0,1.0)

coords1=place_atoms_on_3D_lattice(Npd,bs_2)
coords2=deepcopy(coords1)
atoms=[Atom(σ=1.0,ϵ=1.0,mass=1.0) for i=1:N]
M=[a.mass for a=atoms]
velocities1=init_velocities(T,M,1.0)
velocities2=deepcopy(velocities1)


logger1=SelfDiffusionLogger(coords1)
logger2=SelfDiffusionLogger(coords2)

sys1=System(atoms=atoms,coords=coords1,velocities=velocities1,loggers=Dict(:sd=>logger1),box_size=bs_1,force_units=NoUnits,energy_units=NoUnits)
sys2=System(atoms=atoms,coords=coords2,velocities=velocities2,loggers=Dict(:sd=>logger2),box_size=bs_2,force_units=NoUnits,energy_units=NoUnits)

sim1=LangevinSplitting(T=T,γ=1.0,dt=5e-3,rseed=rseed,splitting="BAOAB")
sim2=deepcopy(sim1)

simulate!(sys1,sim1,2000000)
simulate!(sys2,sim2,2000000)

for i=1:N
    println(sys1.loggers[:sd].self_diffusion_coords[i])
    println(sys2.loggers[:sd].self_diffusion_coords[i])
    println("-----")
end