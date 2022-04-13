include("../molly/MollyExtend.jl")
using .MollyExtend

Npd=10
rseed=2020
N=Npd^3
T=1.0

bs=SVector(3.0,3.0,3.0)

coords1=place_atoms_on_3D_lattice(Npd,bs)
atoms=[Atom(σ=1.0,ϵ=1.0,mass=1.0) for i=1:N]
M=[a.mass for a=atoms]
velocities1=init_velocities(T,M,1.0)


logger1=SelfDiffusionLogger(coords1)


sys1=System(atoms=atoms,coords=coords1,velocities=velocities1,loggers=Dict(:sd=>logger1),box_size=bs,force_units=NoUnits,energy_units=NoUnits)
sim1=LangevinSplitting(T=T,γ=1.0,dt=5e-3,rseed=rseed,splitting="BAOAB")
#sim1=VelocityVerlet(dt=1e-3)
N_steps=10000

W=zeros(N_steps,3)
P=zeros(N_steps,3)
for i=1:N_steps
    simulate!(sys1,sim1,1)
    W[i,:]=first(sys1.coords)
    P[i,:]=first(sys1.loggers[:sd].self_diffusion_coords)
end

using Plots

plot(W[:,1],color=:red)
plot!(P[:,1],color=:red)

plot!(W[:,2],color=:green)
plot!(P[:,2],color=:green)

plot!(W[:,3],color=:blue)
plot!(P[:,3],color=:blue)

savefig("plot.pdf")