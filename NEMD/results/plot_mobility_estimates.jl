using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/NEMD/*.out"
path_end="."
node_color="clustern15"
node_single="clustern16"
run(`scp $node_color:$path_orig $path_end`)
run(`scp $node_single:$path_orig $path_end`)

n_regr=5

file_regex=r"mobility_estimates(.+)_(.+)\.out"

files=readdir()
files=[f for f in files if occursin(file_regex,f)]

ηs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])
Rs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[])


joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
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

for m in ["COLOR","SINGLE"]
  single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  #(m=="COLOR") && (Rs*=1000)
  perm=sortperm(ηs[m])
  ηs[m]=ηs[m][perm]
  Rs[m]=Rs[m][perm]

  a=inv(dot(ηs[m][1:n_regr],ηs[m][1:n_regr]))*dot(ηs[m][1:n_regr],Rs[m][1:n_regr])#least squares slope fit
  println(a)
  scatter!(single_plot,ηs[m],Rs[m],markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(single_plot,x->a*x,0,maximum(ηs[m]),linestyle=:dot,color=:red,label="slope $(round(a,digits=2))")
  scatter!(joint_plot,ηs[m],Rs[m],markershape=:xcross,label=m)
  plot!(joint_plot,x->a*x,0,maximum(ηs[m]),linestyle=:dot,label="slope $(round(a,digits=2))")
  savefig(single_plot,"$(m).pdf")
end
savefig(joint_plot,"joint.pdf")
