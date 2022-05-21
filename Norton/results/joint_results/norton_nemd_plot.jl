using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_nemd="/libre/blasseln/MolecularDynamicsJulia/NEMD/results/"
path_norton="/libre/blasseln/MolecularDynamicsJulia/Norton/results/"

n_linear_regime=10

nemd_regex=r"mobility_estimates(.+)_(.+)\.out"
norton_regex=r"norton_mobility_estimates(.+)_(.+)\.out"

files_nemd=[f for f in readdir(path_nemd) if occursin(nemd_regex,f)]
files_norton=[f for f in readdir(path_norton) if occursin(norton_regex,f)]

ηs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
Rs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])

vs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
dΛs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])

methods=["SINGLE","COLOR","TWO"]
joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
joint_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)

for (i,f)=enumerate(files_nemd)
  println("$i/$(length(files_nemd))")
  flush(stdout)
  (η,method)=match(nemd_regex,f)
  η=parse(Float64,η)
  push!(ηs[method],η)
  file_handle=open(path_nemd*f,"r")
  n_samps=0
  sum_R=0.0
  while !eof(file_handle)
    sum_R+=read(file_handle,Float64)
    n_samps+=1
  end
  close(file_handle)
  push!(Rs[method],sum_R/n_samps)
end

for (i,f)=enumerate(files_norton)
    println("$i/$(length(files_norton))")
    flush(stdout)
    (v,method)=match(norton_regex,f)
    v=parse(Float64,v)
    push!(vs[method],v)
    file_handle=open(path_norton*f,"r")
    n_samps=0
    sum_lambda=0.0
    while !eof(file_handle)
      sum_lambda+=read(file_handle,Float64)
      n_samps+=1
    end
    close(file_handle)
    push!(dΛs[method],sum_lambda/n_samps)
  end

for m in methods
  single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  single_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  perm_nemd=sortperm(ηs[m])

  ηs[m]=ηs[m][perm_nemd]
  Rs[m]=Rs[m][perm_nemd]

  scatter!(single_plot,ηs[m],Rs[m],markershape=:xcross,label="$(m)_T",legend=:topleft)
  scatter!(single_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T",legend=:topleft)
  scatter!(joint_plot,ηs[m],Rs[m],markershape=:xcross,label="$(m)_T")
  scatter!(joint_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T")

  perm_norton=sortperm(dΛs[m])
  dΛs[m]=dΛs[m][perm_norton]
  vs[m]=vs[m][perm_norton]
  
  scatter!(single_plot,dΛs[m],vs[m],markershape=:xcross,label="$(m)_N",legend=:topleft)
  scatter!(single_plot_linear_regime,dΛs[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_N",legend=:topleft)
  scatter!(joint_plot,dΛs[m],vs[m],markershape=:xcross,label="$(m)_N")
  scatter!(joint_plot_linear_regime,dΛs[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_N")

  savefig(single_plot,"$(m).pdf")
  savefig(single_plot_linear_regime,"$(m)_linear_regime.pdf")
end

savefig(joint_plot,"joint_plot.pdf")
savefig(joint_plot_linear_regime,"joint_plot_linear_regime.pdf")