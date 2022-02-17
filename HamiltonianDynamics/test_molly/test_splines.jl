using Molly

include("../molly/custom_cutoffs.jl")

r_c=2.0
r_s=1.5

cutoffs=[NoCutoff(),
DistanceCutoff(r_c),
ShiftedPotentialCutoff(r_c),
ShiftedForceCutoff_(r_c),
CubicSplineCutoff(r_s,r_c)]


r2=1.7^2

R=0.9:0.001:(r_c+0.2)

coord_i=SVector(0.0,0.0,0.0)
atom=Atom(mass=1.0,σ=1.0,ϵ=1.0)
box_size=SVector(10.0,10.0,10.0)

graphs=[Float64[] for i in 1:5]
for r in R
    for (i,cut) in enumerate(cutoffs)
        coord_j=SVector(r,0.0,0.0)
        inter=LennardJones(cutoff=cut,force_units=NoUnits,energy_units=NoUnits)
        push!(graphs[i],Molly.potential_energy(inter,coord_i,coord_j,atom,atom,box_size))
    end
    
end

graphs=hcat(graphs...)



using Plots
pl=plot(xlabel="r/σ",ylabel="U/ϵ",ylims=(-1.2,0.5),dpi=300)
for (i,l) in enumerate(["none","sharp","shift","linear correction","spline"])
    plot!(R,graphs[:,i], label=l)
end
plot!([r_c,r_s],series_type=:vline,color=:black,linestyle=:dot,label="")
vline!([r_c,r_s],color=:black,linestyle=:dash,label="",linewidth=0.5)
savefig(pl,"cutoffs.png")