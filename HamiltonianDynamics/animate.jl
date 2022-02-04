using Plots

function animate_trajectories(sys,filename)

l,l,l=sys.box_size
P=sys.loggers[:position].coords
n_steps=length(P)
N=length(first(P))

anim=@animate for i=1:n_steps
    x,y,z=P[i][1]
    X=[P[i][j][1] for j=2:N]
    Y=[P[i][j][2] for j=2:N]
    Z=[P[i][j][3] for j=2:N]
    plot([x],[y],[z],seriestype=:scatter,color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l))
    plot!(X,Y,Z,seriestype=:scatter,color=:grey,xlims=(0,l),ylims=(0,l),zlims=(0,l))
end

gif(anim,filename,fps=30)

end