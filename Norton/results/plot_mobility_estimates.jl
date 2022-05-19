using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/Norton/*.out"
path_end="."
node_color="clustern23"
node_single="clustern19"
run(`scp $node_color:$path_orig $path_end`)
run(`scp $node_single:$path_orig $path_end`)
v_dict=Dict("COLOR"=>(0.1:0.1:1.0),"SINGLE"=>(0.1:0.1:1.0))

methods=["SINGLE","COLOR"]
joint_plot=plot(ylims=(0,1.0),xlabel="Forcing",ylim="Response",legend=:topleft)
for m in methods
    normal_plot=plot(ylims=(0,1.0),xlabel="Forcing",ylim="Response",legend=:topleft)
    println(m)
       d_lambdas=[]
       vs=v_dict[m]
       for v in vs
        println("\t",v)
        f=open("norton_mobility_estimates$(v)_$(m).out","r")
        n_samps=0
        sum_lambda=0.0
        while !eof(f)
          sum_lambda+=read(f,Float64)
          n_samps+=1
        end
        close(f)
        push!(d_lambdas,sum_lambda/n_samps)
       end
  #(m=="COLOR") && (vs*=1000)
  a=inv(dot(d_lambdas,d_lambdas))*dot(d_lambdas,vs)#least squares slope fit
  println(a)
  scatter(normal_plot,d_lambdas,vs,markershape=:xcross,label=m,color=:blue)
  plot!(normal_plot,x->a*x,0,maximum(d_lambdas),linestyle=:dot,color=:red,label="slope: $(round(a,digits=2))")
  scatter!(joint_plot,d_lambdas,vs,markershape=:xcross,label=m)
  plot!(joint_plot,x->a*x,0,maximum(d_lambdas),linestyle=:dot,label="slope: $(round(a,digits=2))")
  savefig(normal_plot,"$(m).pdf")
end
savefig(joint_plot,"joint_plot.pdf")