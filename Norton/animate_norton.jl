using LinearAlgebra, Plots

include("/home/noeblassel/Documents/stage_CERMICS_2022/molly/MollyExtend.jl")

using .MollyExtend

Npd=6

ρ = 0.7
T = 1.0

r_a = 2.5
r_c = 4.0

Npd = 10
N = Npd^3
γ=1.0
dt=1e-3

L = (N / ρ)^(1 // 3)

box_size = SVector(L, L, L)

v=5.0

coords = place_atoms_on_3D_lattice(Npd,box_size)
atoms = [Atom(σ = 1.0, ϵ = 1.0, mass = 1.0) for i in 1:N]
velocities = init_velocities(T,[a.mass for a=atoms],1.0)

nf = nothing

if 3r_c<L
    global nf = CellListMapNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c, unit_cell = box_size)
else
    global nf = TreeNeighborFinder(nb_matrix = trues(N, N), dist_cutoff = r_c)
end


ff=TwoDriftNEMD(N,1.0).force_field

R(s::System,neighbors=nothing)=v-dot(ff,accelerations(s,neighbors))

inter = LennardJones(cutoff = ShiftedForceCutoff(r_c), nl_only = true, force_units = NoUnits, energy_units = NoUnits)
simulator=NortonTwoDriftSplitting(N=N,dt = dt, γ = γ, T = T, v=v,splitting="BAOAB")
loggers = Dict(:coords => CoordinateLogger(Float64, 1),:lambdas=>GeneralObservableLogger(R,1))

n_eq_steps=5000
n_steps=5000

sys = System(atoms = atoms, coords = coords, velocities = velocities, pairwise_inters = (inter,), box_size = box_size, neighbor_finder = nf, force_units = NoUnits, energy_units = NoUnits, loggers = loggers)

simulate!(sys,simulator,n_eq_steps)
empty!(sys.loggers[:coords].coords)
empty!(sys.loggers[:lambdas].history)
println("equilibriated")
ts=(1:n_steps)*dt
simulate!(sys,simulator,n_steps)
dΛ_hist=sys.loggers[:lambdas].history
plot(ts,dΛ_hist)
savefig("norton_fluctuations_color_drift_exag_.pdf")


println(sys.velocities[1],sys.velocities[2])
println(simulator.F[1],simulator.F[2])

function animate_color_drift(sys,filename)
    l,l,l=sys.box_size
    P=sys.loggers[:coords].coords
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
        
        canvas=plot(xlims=(0,l),ylims=(0,l),zlims=(0,l),showaxis=true,ticks=false,camera=(0,0))
        plot!(canvas,X_red,Y_red,Z_red,seriestype=:scatter,color=:red,label="",markersize=6)
        plot!(canvas,X_blue,Y_blue,Z_blue,seriestype=:scatter,color=:blue,label="",markersize=6)
        plot!(canvas,[0,1],[0,0],[0,0],color=:black,linewidth=3,label="",camera=(0,0))

        drift=plot(xlims=(0,last(ts)),ylims=(minimum(dΛ_hist),maximum(dΛ_hist)),xlabel="t",ylabel="dΛ")
        plot!(drift,ts[1:i],dΛ_hist[1:i],label="")

        plot(canvas,drift,layout=grid(2,1,heights=[0.8,0.2]),size=(900,800))
    end
    mp4(anim,filename,fps=30)

end

function animate_one_drift(sys,filename;ix=rand(1:length(coords[1])))

    l,l,l=sys.box_size
    P=sys.loggers[:coords].coords
    n_steps=length(P)
    N=length(first(P))

    traj_bits_x=[]
    current_traj_bit_x=[]
    traj_bits_y=[]
    current_traj_bit_y=[]
    traj_bits_z=[]
    current_traj_bit_z=[]

    anim=@animate for i=1:n_steps

        canvas=plot(xlims=(0,l),ylims=(0,l),zlims=(0,l),showaxis=true,ticks=false,camera=(0,0))
        drift=plot(xlims=(0,last(ts)),ylims=(minimum(dΛ_hist),maximum(dΛ_hist)),xlabel="t",ylabel="dΛ")

        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end
        x,y,z=popat!(P[i],ix)
        X=[P[i][j][1] for j=1:N-1]
        Y=[P[i][j][2] for j=1:N-1]
        Z=[P[i][j][3] for j=1:N-1]
    
        plot!(canvas,[x],[y],[z],seriestype=:scatter,color=:red,label="",markersize=3)
        plot!(canvas,X,Y,Z,seriestype=:scatter,color=:gray,label="",markersize=2)
        
        for j=1:length(traj_bits_x)
            plot!(canvas,traj_bits_x[j],traj_bits_y[j],traj_bits_z[j],color=:red,label="")
        end
        if i>1
            lx=last(current_traj_bit_x)
            ly=last(current_traj_bit_y)
            lz=last(current_traj_bit_z)

            if norm([x-lx,y-ly,z-lz])>l/2
                if length(current_traj_bit_x)>1
                    plot!(canvas,current_traj_bit_x,current_traj_bit_y,current_traj_bit_z,color=:red,label="")
                end
                
                push!(traj_bits_x,current_traj_bit_x)
                push!(traj_bits_y,current_traj_bit_y)
                push!(traj_bits_z,current_traj_bit_z)
                current_traj_bit_x=[x]
                current_traj_bit_y=[y]
                current_traj_bit_z=[z]
            else
                push!(current_traj_bit_x,x)
                push!(current_traj_bit_y,y)
                push!(current_traj_bit_z,z)
                if length(current_traj_bit_x)>1
                    plot!(canvas,current_traj_bit_x,current_traj_bit_y,current_traj_bit_z,color=:red,label="")
                end        
            end
        else
            current_traj_bit_x=[x]
            current_traj_bit_y=[y]
            current_traj_bit_z=[z]
        end
        plot!(drift,ts[1:i],dΛ_hist[1:i],label="")
        plot(canvas,drift,layout=grid(2,1,heights=[0.8,0.2]),size=(900,800))
    end
    mp4(anim,filename,fps=30)

end
animate_one_drift(sys,"norton_two_drift_exag_1.mp4",ix=1)
animate_one_drift(sys,"norton_two_drift_exag_2.mp4",ix=1)