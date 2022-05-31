include("../molly/MollyExtend.jl")

using Plots, .MollyExtend

Npd=10

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3
γ=1.0
dt=5e-3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

coords = place_atoms_on_3D_lattice(Npd,box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = init_velocities(T,[a.mass for a=atoms],1.0)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
ff=PiecewiseLinearForceProfile(ξ=10.0,L=L)
simulator=LangevinSplitting(dt = dt, γ = γ, T = T,splitting="BAOAB")
loggers = Dict(:coords => CoordinateLogger(Float64, 1))

n_eq_steps=5000
n_steps=5000

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), general_inters=(ff,),box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

simulate!(sys,simulator,n_eq_steps)
empty!(sys.loggers[:coords].coords)

println("equilibriated")

simulate!(sys,simulator,n_steps)


function animate_system(sys,filename)
    l,l,l=sys.box_size
    map_color(y)=RGB(1-y/l,0.0,y/l)#color particles based on y coordinate
    map_radius(z)=4*8/(8+z)
    P=sys.loggers[:coords].coords
    n_steps=length(P)
    N=length(first(P))
    anim=@animate for i=1:n_steps
        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end
        X=[P[i][j][1] for j=1:N]
        Y=[P[i][j][2] for j=1:N]
        Z=[P[i][j][3] for j=1:N]
        
        plot(X,Y,Z,seriestype=:scatter,color=map_color.(Y),label="",xlims=(0,l),ylims=(0,l),zlims=(0,l),showaxis=true,ticks=false,camera=(0,90),msw=0,markersize=map_radius.(Z))
    end
    mp4(anim,filename,fps=30)

end
animate_system(sys,"pwlinear.mp4")