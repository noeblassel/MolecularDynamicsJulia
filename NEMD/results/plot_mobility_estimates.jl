using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/NEMD/*.out"
path_end="."
node_color="clustern15"
node_single="clustern16"
run(`scp $node_color:$path_orig $path_end`)
run(`scp $node_single:$path_orig $path_end`)

η_dict=Dict("COLOR"=>(0.1:0.1:1.0),"SINGLE"=>(0.1:0.1:1.0))

methods=["SINGLE","COLOR"]
joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
for m in methods
    single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
    println(m)
       Rs=[]
       ηs=η_dict[m]
       for η in ηs
        println("\t",η)
        f=open("mobility_estimates$(η)_$(m).out","r")
        n_samps=0
        sum_R=0.0
        while !eof(f)
          sum_R+=read(f,Float64)
          n_samps+=1
        end
        close(f)
        push!(Rs,sum_R/n_samps)
       end
  #(m=="COLOR") && (Rs*=1000)
  a=inv(dot(ηs,ηs))*dot(ηs,Rs)#least squares slope fit
  println(a)
  scatter!(single_plot,ηs,Rs,markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(single_plot,x->a*x,0,last(ηs),linestyle=:dot,color=:red,label="slope $(round(a,digits=2))")
  scatter!(joint_plot,ηs,Rs,markershape=:xcross,label=m)
  plot!(joint_plot,x->a*x,0,last(ηs),linestyle=:dot,label="slope $(round(a,digits=2))")
  savefig(single_plot,"$(m).pdf")
end
savefig(joint_plot,"joint.pdf")
