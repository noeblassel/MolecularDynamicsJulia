using Plots, Statistics,LinearAlgebra

include("../molly/MollyExtend.jl")

ρ = 0.2
T = 1.0

r_c = 4.0

Npd = 20
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
velocities = [SVector(1.0, 1.0, 1.0) for i in 1:N]
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

γ = 1.0
eq_steps=600

dt = 5e-3

seed = 1234

simulator = LangevinTest(dt = dt, γ = γ, T = T, rseed = seed)

loggers = Dict(:velocity => VelocityLogger(Float64,1))

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (), box_size = box_size, force_units = NoUnits, energy_units = NoUnits,loggers=loggers)

simulate!(sys,simulator,eq_steps)
f_chi(x)=sqrt(2/π)*x^2*exp(-x^2/2)

Xrange=0:0.01:5
f_chis=f_chi.(Xrange)

save_every=50
anim=@animate for (i,vs) in enumerate(sys.loggers[:velocity].velocities)
    norms=norm.(vs)
    println(i,"/",eq_steps)
    histogram(norms,normalize=true,xlims=(0,5),ylims=(0,1),dpi=300,label="",bins=100,color=:blue,linecolor=:blue)
    plot!(Xrange,f_chis,label="Maxwell's law",xlabel="velocity",ylabel="likelihood",xlims=(0,5),ylims=(0,1),dpi=300,color=:red)
    (i<=10 || i%save_every==0) && savefig("vel_relaxation/velocities_$(i).png")

end

mp4(anim,"velocity_relaxations.mp4",fps=10)

