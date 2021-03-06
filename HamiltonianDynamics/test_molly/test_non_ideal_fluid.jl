using Molly
using Plots
using ProgressMeter
using Statistics,StaticArrays
using Unitful, UnitfulRecipes


include("../../molly/custom_loggers.jl")
include("../../molly/io.jl")
include("../../molly/custom_observables.jl")
include("../../molly/custom_cutoffs.jl")
include("../../molly/sim_nve_lj.jl")

include("../../utils/place_atoms.jl")
include("../../utils/reduced_units.jl")

Npd=12

Ts=[]
Ps=[]

ρ=0.2
T_range=0.1:0.01:1.0
lrc=0
V=0
sys=nothing
for T in T_range
    sys=sim_lennard_jones_fluid(Npd,ρ,T,1e-3,10000,VelocityVerlet,[],5.0)
    sys.loggers=Dict(:pressure=>PressureLoggerReduced(Float64,1),:temperature=>TemperatureLoggerReduced(Float64,1))
    simulate!(sys,VelocityVerlet(dt=0.1e-3),10000)
    push!(Ts,mean(sys.loggers[:temperature].temperatures))
    push!(Ps,mean(sys.loggers[:pressure].pressures))
    
    if T==first(T_range)
        global lrc=long_range_virial_correction(sys,sys.general_inters[1])
        global V=sys.box_size[1]*sys.box_size[2]*sys.box_size[3]
    end
end

"""scatter(Ts,Ps,label="",xlabel="T/T*",ylabel="P/P*",title="ρ/ρ*=0.8",dpi=300,markershape=:xcross)
savefig("non_ideal_gas.png")"""

file=open("series_data.out","w")
for t in Ts
    print(file,"$(t) ")
end
println(file)
for p in Ps
    print(file,"$(p) ")
end
println(file)
print(file,lrc/3V)
close(file)
