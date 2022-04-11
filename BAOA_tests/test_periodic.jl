include("../molly/MollyExtend.jl")
using .MollyExtend,Molly
include("periodic_potential.jl")


function update_hist!(bins::Vector{Int64},pt::Float64,lims::Tuple{Float64,Float64})
    N=length(bins)
    (lm,lM)=lims
    
    (pt > lM) && return nothing
    (pt < lm) && return nothing

    bin_ix=Int64(ceil(N*(pt-lm)/(lM-lm)))
    bins[bin_ix]+=1
end

println("usage: dt N N_steps(between samples) N_reps(number of samples) splitting")

dt=parse(Float64,ARGS[1])
N=parse(Int64,ARGS[2])
N_steps=parse(Int64,ARGS[3])
N_reps=parse(Int64,ARGS[4])
splitting=ARGS[5]
N_bins=parse(Int64,ARGS[6])

L=1.0

T=1.0
γ=1.0

p_lims=(-10.0,10.0)
q_lims=(0.0,L)

q_hist=zeros(Int64,N_bins)
p_hist=zeros(Int64,N_bins)

coords=range(0.0,L,N)
coords=[SVector{1,Float64}(x) for x in coords]
velocities=zero(coords)
inters=(SinusoidalPotential(),)
atoms=[Atom(mass=1.0,ϵ=1.0,σ=1.0) for i=1:N]

sys=System(coords=coords,velocities=velocities,atoms=atoms,general_inters=inters,box_size=(L,),energy_units=NoUnits,force_units=NoUnits)
sim=LangevinSplitting(dt=dt,γ=γ,T=T,splitting=splitting)
simulate!(sys,sim,N_steps)# equilibration

for i=1:N_reps-1
    println("$(i)/$(N_reps) steps completed")
    flush(stdout)

    Q=[first(q) for q in sys.coords]
    P=[first(p) for p in sys.velocities]

    update_hist!.((q_hist,),Q,(q_lims,))
    update_hist!.((p_hist,),P,(p_lims,))

    simulate!(sys,sim,N_steps)
end

Q=[first(q) for q in sys.coords]
P=[first(p) for p in sys.velocities]

n_samples_q=N_reps*N
n_samples_p=N_reps*N

if isfile("periodic_qs_$(splitting)_$(dt).out")
    r=readlines("periodic_qs_$(splitting)_$(dt).out")
    if length(r)==N_bins+1
        for i=1:N_bins
            q_hist+=parse(Int64,r[i])
        end
        n_samples_q+=parse(Int64,last(r))
    end
end

if isfile("periodic_ps_$(splitting)_$(dt).out")
    r=readlines("periodic_ps_$(splitting)_$(dt).out")
    if length(r)==N_bins+1
        for i=1:N_bins
            p_hist+=parse(Int64,r[i])
        end
        n_samples_p+=parse(Int64,last(r))
    end
end

f=open("periodic_qs_$(splitting)_$(dt).out","w")
print(f,join(q_hist,'\n'))
print(f,"\n$(n_samples_q)")
close(f)

f=open("periodic_ps_$(splitting)_$(dt).out","w")
print(f,join(p_hist,'\n'))
print(f,"\n$(n_samples_p)")
close(f)

