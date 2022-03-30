function animate_trajectories(coords,filename;ix=rand(1:length(coords[1])),camera=(30,30))

    l,l,l=sys.box_size
    P=coords
    n_steps=length(P)
    N=length(first(P))

    traj_bits_x=[]
    current_traj_bit_x=[]
    traj_bits_y=[]
    current_traj_bit_y=[]
    traj_bits_z=[]
    current_traj_bit_z=[]

    anim=@animate for i=1:n_steps
        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end
        x,y,z=popat!(P[i],ix)
        X=[P[i][j][1] for j=1:N-1]
        Y=[P[i][j][2] for j=1:N-1]
        Z=[P[i][j][3] for j=1:N-1]
    
        plot([x],[y],[z],seriestype=:scatter,color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",markersize=3,showaxis=true,ticks=false,camera=camera)
        plot!(X,Y,Z,seriestype=:scatter,color=:gray,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",markersize=2,showaxis=true,ticks=false,camera=camera)
        
        for j=1:length(traj_bits_x)
            plot!(traj_bits_x[j],traj_bits_y[j],traj_bits_z[j],color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",showaxis=true,ticks=false,camera=camera)
        end
        if i>1
            lx=last(current_traj_bit_x)
            ly=last(current_traj_bit_y)
            lz=last(current_traj_bit_z)

            if norm([x-lx,y-ly,z-lz])>l/2
                if length(current_traj_bit_x)>1
                    plot!(current_traj_bit_x,current_traj_bit_y,current_traj_bit_z,color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",showaxis=true,ticks=false,camera=camera)
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
                    plot!(current_traj_bit_x,current_traj_bit_y,current_traj_bit_z,color=:red,xlims=(0,l),ylims=(0,l),zlims=(0,l),label="",showaxis=true,ticks=false,camera=camera)
                end        
            end
        else
            current_traj_bit_x=[x]
            current_traj_bit_y=[y]
            current_traj_bit_z=[z]
        end
    
    end
    mp4(anim,filename,fps=30)

end