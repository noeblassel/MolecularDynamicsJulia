using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/Norton/*.out"
path_end="."
node_color="clustern23"
node_single="clustern19"
run(`scp $node_color:$path_orig $path_end`)
run(`scp $node_single:$path_orig $path_end`)
file_regex=r"norton_mobility_estimates(.+)_(.+)\.out"
files=readdir()
files=[f for f in files if occursin(file_regex,f)]

vs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])
Lambdas=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])
methods=["SINGLE","COLOR"]
joint_plot=plot(ylims=(0,1.0),xlabel="Forcing",ylabel="Response",legend=:topleft)

for (i,f)=enumerate(files)
  println("$i/$(length(files))")
  flush(stdout)
  (v,method)=match(file_regex,f)
  v=parse(Float64,η)
  push!(vs[method],v)
  file_handle=open(f,"r")
  n_samps=0
  sum_lambda=0.0
  while !eof(file_handle)
    sum_lambda+=read(file_handle,Float64)
    n_samps+=1
  end
  close(file_handle)
  push!(Lambdas[method],sum_lambda/n_samps)
end

for m in ["COLOR","SINGLE"]
  single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  #(m=="COLOR") && (Rs*=1000)
  a=inv(dot(Lambdas[m],Lambdas[m]))*dot(Lambdas[m],vs[m])#least squares slope fit
  println(a)
  scatter!(single_plot,Lambdas[m],vs[m],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(single_plot,x->a*x,0,maximum(Lambdas[m]),linestyle=:dot,color=:red,label="slope $(round(a,digits=2))")
  scatter!(joint_plot,Lambdas[m],vs[m],markershape=:xcross,label=m)
  plot!(joint_plot,x->a*x,0,maximum(Lambdas[m]),linestyle=:dot,label="slope $(round(a,digits=2))")
  savefig(single_plot,"$(m).pdf")
end
savefig(joint_plot,"joint.pdf")
savefig(joint_plot,"joint_plot.pdf")