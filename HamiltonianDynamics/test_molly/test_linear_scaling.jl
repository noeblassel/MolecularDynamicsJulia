using Pkg

Pkg.add("Molly")
Pkg.add("ProgressMeter")
Pkg.add("Plots")


using Molly, ProgressMeter, Plots

include("../molly/SimulateLennardJones.jl")
include("../molly/custom_simulators.jl")
include("../molly/custom_loggers.jl")

times_no_nf = []
times_nf_tree = []
times_nf_dist = []
times_nf_cell=[]

ρ = 0.2
T = 0.5
r_c = 3.0
initial_configs = []

Npd_range = 4:12
N_range = Npd_range .^ 3

inter_nl=(LennardJones(cutoff = ShiftedForceCutoff(r_c),nl_only=true,force_units = NoUnits, energy_units = NoUnits),)
inter_no_nl=(LennardJones(cutoff = ShiftedForceCutoff(r_c),force_units = NoUnits, energy_units = NoUnits),)

for Npd in Npd_range
    conf = sim_lennard_jones_fluid(Npd, ρ, T, 5e-3, 10000, VelocityVerlet, [], r_c)
    t_matrix=trues(length(conf), length(conf))

    coords, vels= deepcopy.([conf.coords, conf.velocities])
    sys_no_nf = System(atoms = conf.atoms, general_inters = inter_no_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = NoNeighborFinder())
    push!(times_no_nf, @elapsed sim_lennard_jones_fluid!(sys_no_nf, 5e-3, 10000, VelocityVerlet))

    coords, vels= deepcopy.([conf.coords, conf.velocities])
    nf_tree = TreeNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
    sys_nf_tree = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_tree)
    push!(times_nf_tree, @elapsed sim_lennard_jones_fluid!(sys_nf_tree, 5e-3, 10000, VelocityVerlet))

    coords, vels= deepcopy.([conf.coords, conf.velocities])
    nf_dist = DistanceNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
    sys_nf_dist = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_dist)
    push!(times_nf_dist, @elapsed sim_lennard_jones_fluid!(sys_nf_dist, 5e-3, 10000, VelocityVerlet))

    coords, vels= deepcopy.([conf.coords, conf.velocities])
    nf_cell = CellListMapNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
    sys_nf_cell = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_cell)
    push!(times_nf_cell, @elapsed sim_lennard_jones_fluid!(sys_nf_cell, 5e-3, 10000, VelocityVerlet))
end

plot(N_range, times_no_nf, xlabel = "N", ylabel = "runtime", label = "none", title = "Neighbor finding strategies",xaxis:log,yaxis:log)
plot!(N_range, times_nf_tree, label = "tree")
plot!(N_range,times_nf_dist,label="distance")
plot!(N_range,times_nf_cell,label="cell")

savefig("test_linear_scaling.jl")