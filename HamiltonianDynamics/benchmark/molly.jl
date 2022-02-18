using Pkg

Pkg.add("Molly")
Pkg.add("ProgressMeter")
Pkg.add("Plots")


using Molly, ProgressMeter, Plots

include("../molly/SimulateLennardJones.jl")
include("../molly/custom_simulators.jl")
include("../molly/custom_loggers.jl")
