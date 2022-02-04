using Plots

function animate_trajectories(sys,filename)

P=sys.loggers[:position].coords
n_steps=length(P)
N=length(first(P))

anim=@animate for i=1:n_steps
    x,y,z=P[i][1]
    X=[P[i][j][1] for j=2:N]
    Y=[P[i][j][2] for j=2:N]
    Z=[P[i][j][3] for j=2:N]
    plot([x],[y],[z],seriestype=:scatter,color=:red)
    plot!(X,Y,Z,seriestype=:scatter,color=:grey)
end

gif(anim,filename,fps=30)

end