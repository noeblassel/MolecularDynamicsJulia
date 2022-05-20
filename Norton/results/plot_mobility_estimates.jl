using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/Norton/*.out"
path_end="."
node_color="clustern23"
node_single="clustern19"
run(`scp $node_color:$path_orig $path_end`)
run(`scp $node_single:$path_orig $path_end`)

n_regr=5
n_linear_regime=20

file_regex=r"norton_mobility_estimates(.+)_(.+)\.out"
files=readdir()
files=[f for f in files if occursin(file_regex,f)]

vs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])
Lambdas=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])
methods=["SINGLE","COLOR"]
joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
joint_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)

for (i,f)=enumerate(files)
  println("$i/$(length(files))")
  flush(stdout)
  (v,method)=match(file_regex,f)
  v=parse(Float64,v)
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
  single_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  perm=sortperm(Lambdas[m])

  Lambdas[m]=Lambdas[m][perm]
  vs[m]=vs[m][perm]

  a=inv(dot(Lambdas[m][1:n_regr],Lambdas[m][1:n_regr]))*dot(Lambdas[m][1:n_regr],vs[m][1:n_regr])#least squares slope fit
  println(a)
  scatter!(single_plot,Lambdas[m],vs[m],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  scatter!(single_plot_linear_regime,Lambdas[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(single_plot_linear_regime,x->a*x,linestyle=:dot,color=:red,label="slope $(round(a,digits=2))")
  scatter!(joint_plot,Lambdas[m],vs[m],markershape=:xcross,label=m)
  scatter!(joint_plot_linear_regime,Lambdas[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label=m)
  plot!(joint_plot_linear_regime,x->a*x,linestyle=:dot,label="slope $(round(a,digits=2))")
  savefig(single_plot,"$(m).pdf")
end
savefig(joint_plot,"joint.pdf")
savefig(joint_plot,"joint_plot.pdf")