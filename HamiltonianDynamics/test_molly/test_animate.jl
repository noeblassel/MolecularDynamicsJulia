include("../molly/io.jl")
include("../utils/animate.jl")

using Molly
using ProgressMeter

sys=read_reduced_lj_state("starting_states/T(0.1).out")
sys.loggers=Dict(:coords=>CoordinateLogger(Float64,1))
simulate!(sys,VelocityVerlet(dt=0.005),5000)

animate_trajectories(sys.loggers[:coords].coords,"anim_cold.mp4")