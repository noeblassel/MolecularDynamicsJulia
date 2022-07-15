using Molly, LinearAlgebra, Plots, Random

include("../utils/place_atoms.jl")
include("norton_integrators.jl")

struct NortonSVIntegrator{TF,TG,TDG}
    dt::Float64
    η::Float64
    T::Float64
    γ::Float64
    F::TF # forcing profile
    G::TG # response profile
    dG::TDG # derivative of response profile
end

NortonSVIntegrator(dt::Float64, η::Float64, T::Float64, γ::Float64, F::Function, G::Function, dG::Function) = NortonSVIntegrator{typeof(F),typeof(G),typeof(dG)}(dt, η, T, γ, F, G, dG)

function simulate_norton!(sys::System{D}, sim::NortonSVIntegrator, n_steps; parallel::Bool=true, rng=Random.GLOBAL_RNG) where {D}
    λ_ham_hist = Float64[]
    λ_fd_hist = Float64[]

    N = length(sys)
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt(1 - α^2)

    println(α)
    println(σ)

    noise = zero(sys.velocities)

    accels = accelerations(sys, neighbors; parallel=parallel)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)
    accels_array = reinterpret(reshape, Float64, accels)
    noise_array = reinterpret(reshape, Float64, noise)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    v_y = view(velocities_array, 2, :)
    v_x_perp=view(velocities_array, 2:D,:)
    q_y = view(coords_array, 2, :)
    a_x = view(accels_array, 1, :)
    w_x = view(noise_array, 1, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)
    dG_y = sim.dG.(q_y)

    #compute useful dot products
    FdotG = dot(F_y, G_y)
    GdotA = dot(G_y, a_x)
    GdotG = dot(G_y, G_y)
    FdotF = dot(F_y, F_y)

    # p-wise reprojection to initialize state
    v_x_base = v_x * sim.η / dot(G_y, v_x) #find particular solution on the target manifold
    v_x .-= (dot(G_y, v_x - v_x_base) / dot(G_y, G_y)) * G_y #orthogonal projection on target manifold with respect to the p coordinate

    run_loggers!(sys, neighbors, 0; parallel=parallel)

    for step_n = 1:n_steps
        #B step
        λ_13 = -(GdotA + dot(dG_y, v_x .* v_y)) / FdotG
        sys.velocities .+= accels * sim.dt / 2
        v_x .+= λ_13 * F_y
        v_x_base .= v_x * sim.η / dot(G_y, v_x)
        v_x .-= (dot(G_y, v_x - v_x_base) / GdotG) * G_y

        #A step
        sys.coords .+= sys.velocities * sim.dt
        sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

        accels .= accelerations(sys, neighbors; parallel=parallel)

        F_y .= sim.F.(q_y)
        G_y .= sim.G.(q_y)
        dG_y .= sim.dG.(q_y)

        FdotG = dot(F_y, G_y)
        (isnan(FdotG)) && return λ_ham_hist,λ_fd_hist
        GdotA = dot(G_y, a_x)
        GdotG = dot(G_y, G_y)
        FdotF = dot(F_y, F_y)

        #B step
        λ_23 = -(GdotA + dot(dG_y, v_x .* v_y)) / FdotG
        sys.velocities .+= accels * sim.dt / 2
        v_x .+= λ_23 * F_y
        v_x_base .= v_x * sim.η / dot(G_y, v_x)
        v_x .-= (dot(G_y, v_x - v_x_base) / GdotG) * G_y

        if step_n % 500 == 0 #verify that the system has not NaNed out, else check that response is constant
            println(step_n, " : ", dot(G_y, v_x))
        end
        neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
        run_loggers!(sys, neighbors, step_n; parallel=parallel)

        push!(λ_ham_hist, λ_13)
        push!(λ_ham_hist, λ_23)

        #O_step

         ## dissipation
         P = I - (F_y*G_y')/FdotG
         Pcomp=I-P
        v_x_perp .*= α
        v_x .= (P + α*Pcomp)*v_x

        ## computation of fluctuation term
        noise .= σ * random_velocities(sys, sim.T; rng=rng) # equilibrium fluctuation
        Σ=sqrt(P*P')
        noise .= Σ*noise

        """w_x .-= dot(w_x, F_y) * F_y / FdotF # fluctuation is zero on subspace spanned by F
        w_x .+= (1 - sqrt(FdotF * GdotG / FdotG^2)) * dot(w_x, G_y) * G_y / GdotG # fluctuation is scaled by 1/|cos(θ)| on the subspace spanned by G in F's orthogonal complement (θ is the angle between F and G)
        ## add fluctuation"""
        sys.velocities .+= noise
        v_x_base .= v_x * sim.η / dot(G_y, v_x)
        v_x .-= (dot(G_y, v_x - v_x_base) / GdotG) * G_y
        
        push!(λ_fd_hist,sim.γ * sim.η/FdotG)
        println(dot(G_y,v_x))
    end
    return λ_ham_hist,λ_fd_hist
end

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 5
N = Npd^3
γ = 1.0
dt = 1e-3

L = (N / ρ)^(1 // 3)

box_size = CubicBoundary(L, L, L)

coords = place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

nf = nothing

if 3r_c < L
    global nf = CellListMapNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c, unit_cell=box_size)
else
    global nf = TreeNeighborFinder(nb_matrix=trues(N, N), dist_cutoff=r_c)
end

inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), nl_only=true, force_units=NoUnits, energy_units=NoUnits)
F(y) = sin(2π * y / L)
F_const(y) = (y < L / 2) ? 1.0 : -1.0
F_lin(y) = (y < L / 2) ? 4 * (y - L / 4) / L : -4 * (3L / 4 - y) / L
F_lin_shifted(y) = F_lin((L / 2 -y) % L)

G(y) = sin(2π * y / L) / N
dG(y) = 2π * cos(2π * y / L) / (N * L)

η = 0.05

simulator_eq=LangevinSplitting(dt=dt,temperature=T,friction=1.0,splitting="BAOAB",remove_CM_motion=false)

sys = System(atoms=atoms, coords=coords, velocities=velocities, pairwise_inters=(inter,), boundary=box_size, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, k=1.0)
simulate!(sys,simulator_eq,5000)

simulator = NortonSVIntegrator(dt, η, T, γ, F, G, dG)
coords_obs(sys, args...; kwargs...) = copy(sys.coords)

loggers = (pot=PotentialEnergyLogger(Float64, 1),temp=TemperatureLogger(Float64,1))#

n_eq_steps = 10000

sys = System(atoms=sys.atoms, coords=sys.coords, velocities=sys.velocities, pairwise_inters=(inter,), boundary=sys.boundary, neighbor_finder=nf, force_units=NoUnits, energy_units=NoUnits, loggers=loggers, k=1.0)

println(first(sys.coords))
hist_ham,hist_fd = simulate_norton!(sys, simulator, n_eq_steps)
println(first(sys.coords))
#println(hist_ham," ",hist_fd)

plot(values(sys.loggers.pot)/N,ylims=(-10.0,10.0))
savefig("pot_const.pdf")
plot(values(sys.loggers.temp),ylims=(0.0,10.0))
savefig("temp_const.pdf")

function animate_system(sys, filename)
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
        Z = [P[i][j][3] for j = 1:N]

        plot(X, Y, seriestype=:scatter, color=map_color.(F.(Y)), label="", xlims=(0, l), ylims=(0, l), zlims=(0, l), showaxis=true, ticks=false, camera=(0, 90), msw=0, markersize=3)
    end
    mp4(anim, filename, fps=30)
end

#animate_system(sys,"const_full.mp4")

