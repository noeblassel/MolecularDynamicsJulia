#!/libre/blasseln/julia-1.7.2/bin/julia

using Pkg

Pkg.add("Molly")
Pkg.add("Plots")
Pkg.add("ProgressMeter")
Pkg.add("Statistics")
Pkg.add("StaticArrays")

using Molly
using Plots
using ProgressMeter
using Statistics,StaticArrays

include("../molly/custom_loggers.jl")
include("../molly/io.jl")
include("../utils/animate.jl")

N=length(sys)
ρ=0.01
Ts=[]
Ps=[]
for T in 0.01:0.1:10.0
    sys=read_reduced_lj_state("starting_states/T($(T)).out")
    sys.loggers=Dict(:pressure=>PressureLoggerReduced(Float64,1),:temperature=>TemperatureLogger(Float64,1))
    simulate!(sys,VelocityVerlet(dt=0.005),5000)
    push!(Ts,mean(sys.loggers[:temperature].temperatures))
    push!(Ps,mean(sys.loggers[:pressure].pressures))
end
    
scatter(Ts,Ps)
T_range=0.01:0.1:10.0
plot!(T_range,ρ*T_range,label="Ideal gas law",xlabel="T",ylabel="P",title="Ideal gas regime, ρ=0.01")
savefig("ideal_gas.png")
