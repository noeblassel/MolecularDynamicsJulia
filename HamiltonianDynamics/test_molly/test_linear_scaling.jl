using Pkg

Pkg.add("Molly")
Pkg.add("ProgressMeter")
Pkg.add("Plots")


using Molly, ProgressMeter, Plots

include("../molly/SimulateLennardJones.jl")
include("../molly/custom_simulators.jl")
include("../molly/custom_loggers.jl")
include("../molly/custom_cutoffs.jl")

times_no_nf = []
times_nf_tree = []
times_nf_dist = []
times_nf_cell=[]

ρ = 0.2
T = 1.0
r_c = 3.0

Npd_range = 4:15
N_range = Npd_range .^ 3

n_samps=10

inter_nl=(LennardJones(cutoff = ShiftedForceCutoff_fixed(r_c),nl_only=true,force_units = NoUnits, energy_units = NoUnits),)
inter_no_nl=(LennardJones(cutoff = ShiftedForceCutoff_fixed(r_c),force_units = NoUnits, energy_units = NoUnits),)

for Npd in Npd_range
    conf = sim_lennard_jones_fluid(Npd, ρ, T, 5e-3, 10000, VelocityVerlet, [], r_c)
    t_matrix=trues(length(conf), length(conf))

    t_no_nf=0
    t_nf_tree=0
    t_nf_dist=0
    t_nf_cell=0


    for i=1:n_samps
        coords, vels= deepcopy.([conf.coords, conf.velocities])
        sys_no_nf = System(atoms = conf.atoms, general_inters = inter_no_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = NoNeighborFinder())
        t_no_nf+=@elapsed sim_lennard_jones_fluid!(sys_no_nf, 5e-3, 10000, VelocityVerlet)

        coords, vels= deepcopy.([conf.coords, conf.velocities])
        nf_tree = TreeNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
        sys_nf_tree = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_tree)
        t_nf_tree+= @elapsed sim_lennard_jones_fluid!(sys_nf_tree, 5e-3, 10000, VelocityVerlet)

        coords, vels= deepcopy.([conf.coords, conf.velocities])
        nf_dist = DistanceNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
        sys_nf_dist = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_dist)
        t_nf_dist+= @elapsed sim_lennard_jones_fluid!(sys_nf_dist, 5e-3, 10000, VelocityVerlet)

        coords, vels= deepcopy.([conf.coords, conf.velocities])
        nf_cell = CellListMapNeighborFinder(nb_matrix = t_matrix, n_steps = 1, dist_cutoff = r_c)
        sys_nf_cell = System(atoms = conf.atoms, general_inters = inter_nl, coords = coords, velocities = vels, box_size = conf.box_size, energy_units = NoUnits, force_units = NoUnits, neighbor_finder = nf_cell)
        t_nf_cell+= @elapsed sim_lennard_jones_fluid!(sys_nf_cell, 5e-3, 10000, VelocityVerlet)
    end

    push!(times_no_nf,t_no_nf/n_samps)
    push!(times_nf_dist,t_nf_dist/n_samps)
    push!(times_nf_tree,t_nf_tree/n_samps)
    push!(times_nf_cell,t_nf_cell/n_samps)
end

f=open("molly_no_nf.txt","w")
for t=times_no_nf
    print(f,"$(t) ")
end
close(f)

f=open("molly_nf_dist.txt","w")
for t=times_nf_dist
    print(f,"$(t) ")
end
close(f)

f=open("molly_nf_tree.txt","w")
for t=times_nf_tree
    print(f,"$(t) ")
end
close(f)

f=open("molly_nf_cell.txt","w")
for t=times_nf_cell
    print(f,"$(t) ")
end
close(f)

f=open("n_range","w")
for n in N_range
print(f,"$(n) ")
end
close(f)
