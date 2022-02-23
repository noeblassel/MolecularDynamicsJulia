include("../../molly/MollyExtend.jl")

sys=read_reduced_lj_state("starting_states/T(0.1).out")
sys.loggers=Dict(:virial=>VirialLogger(Float64,1),:kinetic_energy=>KineticEnergyLoggerNoDims(Float64,1),:pressure=>PressureLoggerReduced(Float64,1))
simulate!(sys,VelocityVerlet(dt=0.005),5000)

plot(sys.loggers[:virial].energies/length(sys))
plot!(sys.loggers[:kinetic_energy].energies/length(sys))
savefig("virial_contribution.png")

plot(sys.loggers[:pressure].pressures)
