using Molly
using Plots
using ProgressMeter
using Statistics,StaticArrays
using UnitfulRecipes

include("../molly/custom_loggers.jl")
include("../molly/io.jl")
include("../utils/animate.jl")
include("../utils/ReducedUnits.jl")

N=length(sys)
Ts=[]
Ps=[]
ρ=0.5
T_range=0.01:0.01:0.5
for T in T_range
    sys=read_reduced_lj_state("starting_states/T($(round(T,digits=2))).out")
    sys.loggers=Dict(:pressure=>PressureLoggerReduced(Float64,1),:temperature=>TemperatureLogger(Float64,1))
    simulate!(sys,VelocityVerlet(dt=0.005),10000)
    push!(Ts,mean(sys.loggers[:temperature].temperatures))
    push!(Ps,mean(sys.loggers[:pressure].pressures))
end

Ts=get_physical_temperature.(:Ar,Ts)
Ps=get_physical_pressure.(:Ar,Ps)

scatter(Ts,Ps,label="",xlabel="T",ylabel="P",title="Non-ideal regime, ρ=$(round(typeof(1.0u"mol * m^-3"),get_physical_density(:Ar,ρ)))")
savefig("non_ideal_gas.png")