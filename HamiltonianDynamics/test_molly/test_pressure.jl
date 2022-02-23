using Plots,Statistics

include("../../molly/MollyExtend.jl")


#reasonable values for (T,ρ): (0.98,0.8)

T_ini=0.98
ρ=0.8

sys=sim_lennard_jones_fluid_nve(12,ρ,T_ini,5e-3,20000,VelocityVerlet,[(:pressure,1),(:temperature,1)],4.0)
plot(sys.loggers[:pressure].pressures)
plot!(sys.loggers[:temperature].temperatures)