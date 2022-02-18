using Molly,LinearAlgebra

include("../molly/custom_cutoffs.jl")

r_c=2.5
r_s=1.5

cutoffs=[NoCutoff(),
DistanceCutoff(r_c),
ShiftedPotentialCutoff(r_c),
ShiftedForceCutoff_fixed(r_c), #fixed version
CubicSplineCutoff(r_s,r_c)]


R=0.9:0.001:(r_c+0.2)

coord_i=SVector(0.0,0.0,0.0)
atom=Atom(mass=1.0,σ=1.0,ϵ=1.0)
box_size=SVector(10.0,10.0,10.0)

graphs=[Float64[] for i in 1:5]
for r in R
    for (i,cut) in enumerate(cutoffs)
        coord_j=SVector(r,0.0,0.0)
        inter=LennardJones(cutoff=cut,force_units=NoUnits,energy_units=NoUnits)
        push!(graphs[i],Molly.potential(inter,coord_i,coord_j,atom,atom,box_size))
    end
    
end

graphs=hcat(graphs...)



using Plots
pl=plot(xlabel="r/σ",ylims=(0.0,5.0),ylabel="F/F*",dpi=300)
for (i,l) in enumerate(["none","sharp","shift","linear correction","spline"])
    plot!(R,graphs[:,i], label=l)
end
vline!([r_c,r_s],color=:black,linestyle=:dash,label="",linewidth=0.5)
savefig(pl,"forces_1.5_2.5.png")