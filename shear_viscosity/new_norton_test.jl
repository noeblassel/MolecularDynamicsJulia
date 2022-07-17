using Molly, Plots, Statistics

include("../utils/place_atoms.jl")
include("norton_integrators.jl")

ρ = 0.7
T = 0.8

r_c = 2.5

Npd = 10
N = Npd^3
γ = 1.0
dt = 5e-3

L = (N / ρ)^(1 // 3)

box_size = CubicBoundary(L, L, L)

coords = place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

max_speed=10.0*sqrt(T)
n_steps_neighbors=floor(Int64,0.2*r_c/(dt*max_speed))

nf = nothing

if 3.6r_c < L
    global nf = CellListMapNeighborFinder(nb_matrix=trues(N, N),n_steps=n_steps_neighbors,dist_cutoff=r_c, unit_cell=box_size)
else
    global nf = TreeNeighborFinder(nb_matrix=trues(N, N), n_steps=n_steps_neighbors,dist_cutoff=r_c)
end

inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)
F(y) = sin(2π * y / L)
F_const(y) = (y < L / 2) ? 1.0 : -1.0
F_lin(y) = (y < L / 2) ? 4 * (y - L / 4) / L : 4 * (3L / 4 - y) / L

G(y) = sin(2π * y / L) / N

η = 0.1

#simulator_eq=LangevinSplitting(dt=dt,temperature=T,friction=1.0,splitting="BAOAB",remove_CM_motion=false)

#simulate!(sys,simulator_eq,5000)

simulator = NortonSVIntegrator(dt, η, T, γ, F_const, G)
coords_obs(sys, args...; kwargs...) = copy(sys.coords)

loggers = (temp=TemperatureLogger(Float64,1),)

n_eq_steps = 5000
sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), boundary=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits,loggers=loggers, k=1.0)

println(maximum(abs.(reinterpret(Float64,sys.velocities))))
println(first(sys.coords))
hist_forcing = simulate!(sys, simulator, n_eq_steps)
println(first(sys.coords))
#println(hist_ham," ",hist_fd)

plot(values(sys.loggers.temp))
savefig("temp.pdf")
plot(hist_forcing)
savefig("forcing.pdf")


println("avg forcing: ",mean(hist_forcing))
println(maximum(abs.(reinterpret(Float64,sys.velocities))))

function animate_system(sys, filename,func)
    l, l, l = sys.boundary.side_lengths
    map_color(y) = RGB(1 - (y + 1) / 2, 0.0, (y + 1) / 2)#color particles based on y coordinate
    P = sys.loggers.coords.history
    n_steps = length(P)
    println(n_steps)
    N = length(first(P))
    anim = @animate for i = 1:n_steps
        if i % 100 == 0
            println("Frame $(i)/$(n_steps)")
        end
        X = [P[i][j][1] for j = 1:N]
        Y = [P[i][j][2] for j = 1:N]

        plot(X, Y, seriestype=:scatter, color=map_color.(func.(Y)), label="", xlims=(0, l), ylims=(0, l), showaxis=true, ticks=false, msw=0, markersize=3)
    end
    mp4(anim, filename, fps=30)
end

