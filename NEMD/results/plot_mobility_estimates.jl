using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/NEMD/*.out"
path_end="."

for node in [15,16]
  run(`scp clustern$node:$path_orig $path_end`)
end

n_regr=5
n_linear_regime=19

file_regex=r"mobility_estimates(.+)_(.+)\.out"
files=readdir()
files=[f for f in files if occursin(file_regex,f)]
η s=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
Rs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
methods=["SINGLE","COLOR","TWO"]
joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
joint_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)

for (i,f)=enumerate(files)
  println("$i/$(length(files))")
  flush(stdout)
  (η,method)=match(file_regex,f)
  η=parse(Float64,η)
  push!(ηs[method],η)
  file_handle=open(f,"r")
  n_samps=0
  sum_R=0.0
  while !eof(file_handle)
    sum_R+=read(file_handle,Float64)
    n_samps+=1
  end
  close(file_handle)
  push!(Rs[method],sum_R/n_samps)
end

for m in methods
  single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  single_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  perm=sortperm(ηs[m])

  ηs[m]=ηs[m][perm]
  Rs[m]=Rs[m][perm]

  a=inv(dot(ηs[m][1:n_regr],ηs[m][1:n_regr]))*dot(ηs[m][1:n_regr],Rs[m][1:n_regr])#least squares slope fit
  println(a)
  scatter!(single_plot,ηs[m],Rs[m],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  scatter!(single_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(single_plot_linear_regime,x->a*x,linestyle=:dot,color=:red,label="slope $(round(a,digits=2))")
  scatter!(joint_plot,ηs[m],Rs[m],markershape=:xcross,label=m)
  scatter!(joint_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],markershape=:xcross,label=m)
  plot!(joint_plot_linear_regime,x->a*x,linestyle=:dot,label="slope $(round(a,digits=2))")
  savefig(single_plot,"$(m).pdf")
  savefig(single_plot_linear_regime,"$(m)_linear_regime.pdf")
end

savefig(joint_plot,"joint_plot.pdf")
savefig(joint_plot_linear_regime,"joint_plot_linear_regime.pdf")