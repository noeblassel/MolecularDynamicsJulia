using Molly,ProgressMeter,Plots,Statistics

include("../utils/ReducedUnits.jl")
include("../utils/PlaceAtoms.jl")

include("../molly/custom_cutoffs.jl")
include("../molly/custom_loggers.jl")
include("../molly/custom_observables.jl")
include("../molly/SimulateLennardJones.jl")

#reasonable values for (T,ρ): (0.98,0.8)

T=0.98
ρ=0.8

sys=sim_lennard_jones_fluid(12,ρ,T,5e-3,20000,VelocityVerlet,[(:pressure,1),(:temperature,1)],4.0)
plot(sys.loggers[:pressure].pressures)
plot!(sys.loggers[:temperature].temperatures)