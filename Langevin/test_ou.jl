using Plots, Statistics

include("../molly/MollyExtend.jl")

ρ = 0.2
T = 0.3

r_a = 2.5
r_c = 4.0

Npd = 6
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
velocities = [reduced_velocity_lj(T) for i in 1:N]
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff_fixed(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing

try
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
catch
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end


γ = 1.0
eq_steps=10000
samp_steps = 30000

dt = 5e-3

seed = 1234

simulator = LangevinTest(dt = dt, γ = γ, T = T, rseed = seed)

sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), general_inters = (), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits)

simulate!(sys, simulator, eq_steps)

loggers = Dict(:velocities => VelocityLogger(Float64, 1))

sys = System(atoms = sys.atoms, coords = sys.coords, velocities = sys.velocities, general_inters = (), box_size = sys.box_size, neighbor_finder = sys.neighbor_finder, force_units = NoUnits, energy_units = NoUnits,loggers=loggers)

simulate!(sys,simulator,samp_steps)

t_range = 0:dt:(samp_steps*dt-dt)


q = zeros(samp_steps, 3N)

for t in 1:samp_steps   
    for i in 1:N
        q[t, 3*(i-1)+1:3*(i-1)+3] = sys.loggers[:velocities].velocities[t][i]
    end
end

E_hat = mean(q, dims = 2)
V_hat = var(q, dims = 2)

#ok for mean and variance

