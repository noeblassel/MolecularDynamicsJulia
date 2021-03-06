using Molly
using Plots
using ProgressMeter
using Statistics

include("../../molly/custom_loggers.jl")
include("../../molly/io.jl")
include("../../utils/animate.jl")


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
    
scatter(Ts,Ps,markershape=:xcross)
T_range=0.01:0.1:10.0
plot!(T_range,ρ*T_range,label="Ideal gas law",xlabel="T",ylabel="P",title="Ideal gas regime, ρ=0.01",linestyle=:dot)
savefig("ideal_gas.png")
