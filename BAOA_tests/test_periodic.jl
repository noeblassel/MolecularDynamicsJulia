include("../molly/MollyExtend.jl")
using .MollyExtend,Molly
include("periodic_potential.jl")

println("usage: dt N N_steps(between samples) N_reps(number of samples) splitting")

dt=parse(Float64,ARGS[1])
N=parse(Int64,ARGS[2])
N_steps=parse(Int64,ARGS[3])
N_reps=parse(Int64,ARGS[4])
splitting=ARGS[5]

L=10.0

T=1.0
γ=1.0

coords=range(0.0,L,N)
coords=[SVector{1,Float64}(x) for x in coords]
velocities=zero(coords)
inters=(SinusoidalPotential(),)
atoms=[Atom(mass=1.0,ϵ=1.0,σ=1.0) for i=1:N]

sys=System(coords=coords,velocities=velocities,atoms=atoms,general_inters=inters,box_size=(L,),energy_units=NoUnits,force_units=NoUnits)
sim=LangevinSplitting(dt=2e-2,γ=γ,T=T,splitting=splitting)
simulate!(sys,sim,N_steps)# equilibration

for i=1:N_reps-1

    println("$(i)/$(N_reps) steps completed")
    flush(stdout)

    Q=[first(q) for q in sys.coords]
    P=[first(v) for v in sys.velocities]
    
    f=open("periodic_qs_$(splitting)_$(dt).out","a")
    print(f,join(Q,'\n'))
    print(f,'\n')
    close(f)

    f=open("periodic_ps_$(splitting)_$(dt).out","a")
    print(f,join(P,'\n'))
    print(f,'\n')
    close(f)
    simulate!(sys,sim,N_steps)
end

Q=[first(q) for q in sys.coords]
P=[first(v) for v in sys.velocities]

f=open("periodic_qs_$(splitting)_$(dt).out","a")
print(f,join(Q,'\n'))
print(f,'\n')
close(f)

f=open("periodic_ps_$(splitting)_$(dt).out","a")
print(f,join(P,'\n'))
print(f,'\n')
close(f)



"""Z=12.660658777519801 #numerically computed normalization constant for dq=1e-6 with a trapezoidal rule
ν(q)=exp(-sin(2π*q/L))/Z

Q=[first(q) for q in sys.coords]
using Plots
histogram(Q,normalize=true)
plot!(ν,0,L)
"""






