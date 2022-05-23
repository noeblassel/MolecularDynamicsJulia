using Plots, LinearAlgebra, Statistics

"""Estimates the asymptotic variance for a correlated time-series based on a block-averaging method."""
function asymptotic_var(v::Vector{Float64})
    data=copy(v)
    avg=mean(v)
    data .-= avg
    N_steps=floor(Int64,log2(length(data)))
    data=data[1:2^N_steps]
    L=1
    N=length(data)
    K=N
    max_var=0
        while K>2
            new_var=var(data)*N*inv(K)
            (new_var<max_var) && break
            max_var=new_var
            K >>=1
            data=[0.5*(data[i]+data[i+1]) for i=1:2:K-1]
        end
    return max_var
end

path_nemd="/libre/blasseln/MolecularDynamicsJulia/NEMD/results/"
path_norton="/libre/blasseln/MolecularDynamicsJulia/Norton/results/"

n_linear_regime=10

nemd_regex=r"mobility_estimates(.+)_(.+)\.out"

files_nemd=[f for f in readdir(path_nemd) if occursin(nemd_regex,f)]

ηs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
vars_nemd=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])

methods=["SINGLE","COLOR","TWO"]
joint_plot=plot(xlabel="Forcing",ylabel="Asymptotic variance",legend=:topleft)
joint_plot_linear_regime=plot(xlabel="Forcing",ylabel="Asymptotic variance",legend=:topleft)

for (i,f)=enumerate(files_nemd)
  println("$i/$(length(files_nemd))")
  flush(stdout)
  (η,method)=match(nemd_regex,f)
  η=parse(Float64,η)
  push!(ηs[method],η)
  file_handle=open(path_nemd*f,"r")
  data=Float64[]

  while !eof(file_handle)
    push!(data,read(file_handle,Float64))
  end
  close(file_handle)
  push!(vars_nemd[method],asymptotic_var(data))
end

for m in methods
  single_plot=plot(xlabel="Forcing",ylabel="Asymptotic variance",legend=:topleft)
  single_plot_linear_regime=plot(xlabel="Forcing",ylabel="Asymptotic Variance",legend=:topleft)
  perm_nemd=sortperm(ηs[m])

  ηs[m]=ηs[m][perm_nemd]
  vars_nemd[m]=vars_nemd[m][perm_nemd]

  scatter!(single_plot,ηs[m],vars_nemd[m],markershape=:xcross,label="$(m)_T",markersize=2)
  scatter!(single_plot_linear_regime,ηs[m][1:n_linear_regime],vars_nemd[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T",legend=:topleft,markersize=2)
  scatter!(joint_plot,ηs[m],vars_nemd[m],markershape=:xcross,label="$(m)_T",markersize=2)
  scatter!(joint_plot_linear_regime,ηs[m][1:n_linear_regime],vars_nemd[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T",markersize=2)

  savefig(single_plot,"vars_$(m).pdf")
  savefig(single_plot_linear_regime,"vars_$(m)_linear_regime.pdf")
end

savefig(joint_plot,"vars_joint_plot.pdf")
savefig(joint_plot_linear_regime,"vars_joint_plot_linear_regime.pdf")