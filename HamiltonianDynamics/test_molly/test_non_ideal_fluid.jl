using Pkg

Pkg.add.(split("Molly Plots ProgressMeter Statistics StaticArrays Unitful UnitfulRecipes"))

using Molly
using Plots
using ProgressMeter
using Statistics,StaticArrays
using Unitful, UnitfulRecipes


include("../molly/custom_loggers.jl")
include("../molly/io.jl")
include("../molly/custom_observables.jl")
include("../molly/custom_cutoffs.jl")
include("../molly/SimulateLennardJones.jl")

include("../utils/PlaceAtoms.jl")
include("../utils/ReducedUnits.jl")

Npd=12

Ts=[]
Ps=[]
LRCs=[]

ρ=0.2
T_range=0.2:0.02:2

for T in T_range
    sys=sim_lennard_jones_fluid(Npd,ρ,T,5e-3,5000,VelocityVerlet,[],4.0)
    sys.loggers=Dict(:pressure=>PressureLoggerReduced(Float64,1),:temperature=>TemperatureLoggerReduced(Float64,1),:lrc=>LRCVirialLogger(Float64,1,sys.general_inters[1]))
    simulate!(sys,VelocityVerlet(dt=0.005),10000)
    push!(Ts,mean(sys.loggers[:temperature].temperatures))
    push!(Ps,mean(sys.loggers[:pressure].pressures))
    push!(LRCs,mean(sys.loggers[:lrc].LRCs))
end

V=Npd^3/ρ

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
for c in LRCs
    print(file,"$(c/3V) ")
end
close(file)