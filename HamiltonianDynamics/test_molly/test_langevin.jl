using Molly, ProgressMeter

include("../molly/custom_simulators.jl")
include("../molly/custom_observables.jl")
include("../molly/custom_loggers.jl")
include("../molly/custom_cutoffs.jl")

include("../utils/ReducedUnits.jl")
include("../utils/PlaceAtoms.jl")

ρ = 0.2
T = 0.5

r_a = 2.5
r_c = 4.0

Npd = 6
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
velocities = [SVector(1.0, 1.0, 1.0) for i in 1:N]
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]

inter = LennardJones(cutoff = CubicSplineCutoff(r_a, r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)

nf = nothing

try
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
catch
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

loggers = Dict(:velocities => VelocityLogger(Float64, 1))

sys = System(atoms = atoms, coords = coords, velocities = deepcopy(velocities), general_inters = (), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)


γ = 0.2
n_steps = 10000
dt = 5e-3

seed = 1234

simulator = LangevinTest(dt = dt, γ = γ, T = T, rseed = seed)
@time simulate!(sys, simulator, n_steps)

t_range = 0:dt:(n_steps*dt-dt)
E(t) = exp(-γ * t)
E_t = E.(t_range)
V(t) = T * (1 - E(2t))
V_t = V.(t_range)

using Plots, Statistics

empirical_E = []
empirical_σ2 = []

for i in 1:n_steps
    v = sys.loggers[:velocities].velocities[i]
    print(typeof(v))
    break
    push!(empirical_E, mean(mean(v[i])))
    push!(empirical_σ2, var(v))
end