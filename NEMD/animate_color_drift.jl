using Plots,LinearAlgebra

include("../molly/MollyExtend.jl")

ρ = 0.2
T = 1.25

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_lattice(Npd, box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = [reduced_velocity_lj(T,atoms[i].mass) for i in 1:N]

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
NEMD_forcing=ColorDriftNEMD(N,20.0)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

loggers = Dict(:coords => CoordinateLogger(Float64, 1))

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,),general_inters=(NEMD_forcing,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)


γ = 1.0
n_steps = 2000
dt = 5e-3

simulator = LangevinSplitting(dt = dt, γ = γ, T = T,splitting="BAOAB")
simulate!(sys,simulator,n_steps)


function animate_color_drift(coords,filename)

    l,l,l=sys.box_size
    P=coords
    n_steps=length(P)
    N=length(first(P))
    
    anim=@animate for i=1:n_steps
        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end

        X_red=[P[i][j][1] for j=1:2:N]
        Y_red=[P[i][j][2] for j=1:2:N]
        Z_red=[P[i][j][3] for j=1:2:N]
    
        X_blue=[P[i][j][1] for j=2:2:N]
        Y_blue=[P[i][j][2] for j=2:2:N]
        Z_blue=[P[i][j][3] for j=2:2:N]
        
        plot(X_red,Y_red,Z_red,seriestype=:scatter,color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",markersize=6,showaxis=true,ticks=false,size=(1024,768),camera=(0,0))
        plot!(X_blue,Y_blue,Z_blue,seriestype=:scatter,color=:blue,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",markersize=6,showaxis=true,ticks=false,camera=(0,0))
        plot!([0,1],[0,0],[0,0],color=:black,linewidth=3,label="",camera=(0,0))
    end
    mp4(anim,filename,fps=30)

end
animate_color_drift(sys.loggers[:coords].coords,"color_drift.mp4")