using Plots,LinearAlgebra, Statistics
#TODO add error bars using block averaging
γ=1.0
"""Estimates the asymptotic variance for a correlated time-series based on a block-averaging method."""
function asymptotic_var(v::Vector{Float64})
    data=copy(v)
    avg=mean(v)
    data .-= avg
    N_steps=floor(Int64,log2(length(data)))
    data=data[1:2^N_steps]#crop series to largest possible power of two
    L=1
    N=length(data)
    K=N
    max_var=0.0
        while K>1000
            max_var=max(max_var,varm(data,0.0)*N*inv(K))
            K >>=1
            data=(data[1:2:end]+data[2:2:end])/2
        end
    return max_var
end


path_nemd="/libre/blasseln/MolecularDynamicsJulia/NEMD/results/"
path_norton="/libre/blasseln/MolecularDynamicsJulia/Norton/results/"

n_linear_regime=10

nemd_regex=r"mobility_estimates(.+)_(.+)\.out"
norton_regex=r"norton_mobility_estimates(.+)_(.+)\.out"

files_nemd=[f for f in readdir(path_nemd) if occursin(nemd_regex,f)]
files_norton=[f for f in readdir(path_norton) if occursin(norton_regex,f)]

ηs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
Rs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
asymptotic_vars_nemd=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
error_bars_nemd=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
n_steps_nemd=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])

vs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
dΛs=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
asymptotic_vars_norton=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
error_bars_norton=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])
n_steps_norton=Dict("COLOR"=>Float64[],"SINGLE"=>Float64[],"TWO"=>Float64[])

methods=["SINGLE","COLOR","TWO"]
joint_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
joint_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
plot_asympt_var=plot(xlabel="Forcing",ylabel="Asymptotic_variance",legend=:topright,xaxis=:log,yaxis=:log)
plot_asympt_var_linear_regime=plot(xlabel="Forcing",ylabel="Asymptotic_variance",legend=:topright,xaxis=:log,yaxis=:log)
plot_nsteps=plot(xlabel="Forcing",ylabel="n_steps")

for (i,f)=enumerate(files_nemd)
  println("$i/$(length(files_nemd))")
  flush(stdout)
  (η,method)=match(nemd_regex,f)
  η=parse(Float64,η)
  push!(ηs[method],η)
  file_handle=open(path_nemd*f,"r")
  n_samps=0
  sum_R=0.0
  data_pts=Float64[]

  while !eof(file_handle)
    push!(data_pts,read(file_handle,Float64))
  end
  close(file_handle)
  σ2=asymptotic_var(data_pts)
  push!(Rs[method],mean(data_pts)) #finite difference linear response estimator
  push!(n_steps_nemd[method],length(data_pts))
  push!(error_bars_nemd[method],sqrt(σ2/length(data_pts)))
  push!(asymptotic_vars_nemd[method],σ2/η^2)
end

for (i,f)=enumerate(files_norton)
    println("$i/$(length(files_norton))")
    flush(stdout)
    (v,method)=match(norton_regex,f)
    v=parse(Float64,v)
    push!(vs[method],v)
    file_handle=open(path_norton*f,"r")
    data_pts=Float64[]

    while !eof(file_handle)
      push!(data_pts,read(file_handle,Float64))
    end

    close(file_handle)

    σ2_ergodic_mean=asymptotic_var(data_pts .- γ*v)
    denom=mean(data_pts)
    push!(dΛs[method],denom) #norton linear response estimator
    push!(n_steps_norton[method],length(data_pts))
    σ2=v^2*σ2_ergodic_mean/denom^4 # delta method
    push!(error_bars_norton[method],sqrt(σ2_ergodic_mean/length(data_pts)))
    push!(asymptotic_vars_norton[method],σ2)
  end

for m in methods
  single_plot=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  single_plot_linear_regime=plot(xlabel="Forcing",ylabel="Response",legend=:topleft)
  perm_nemd=sortperm(ηs[m])

  ηs[m]=ηs[m][perm_nemd]
  Rs[m]=Rs[m][perm_nemd]
  n_steps_nemd[m]=n_steps_nemd[m][perm_nemd]
  error_bars_nemd[m]=error_bars_nemd[m][perm_nemd]
  asymptotic_vars_nemd[m]=asymptotic_vars_nemd[m][perm_nemd]

  scatter!(single_plot,ηs[m],Rs[m],markershape=:xcross,label="$(m)_T",yerror=error_bars_nemd[m],markersize=2,msc=:auto)
  scatter!(single_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],yerror=error_bars_nemd[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T",markersize=2,msc=:auto)
  scatter!(joint_plot,ηs[m],Rs[m],markershape=:xcross,yerror=error_bars_nemd[m],msc=:auto,label="$(m)_T",markersize=2)
  scatter!(joint_plot_linear_regime,ηs[m][1:n_linear_regime],Rs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_T",markersize=2,yerror=error_bars_nemd[m][1:n_linear_regime],msc=:auto)
  scatter!(plot_asympt_var,ηs[m],asymptotic_vars_nemd[m],markershape=:xcross,markersize=2,label="$(m)_T")
  scatter!(plot_asympt_var_linear_regime,ηs[m][1:n_linear_regime],asymptotic_vars_nemd[m][1:n_linear_regime],markershape=:xcross,markersize=2,label="$(m)_T")
  scatter!(plot_nsteps,ηs[m],n_steps_nemd[m],markershape=:xcross,markersize=2,label="$(m)_T")

  perm_norton=sortperm(dΛs[m])
  dΛs[m]=dΛs[m][perm_norton]
  vs[m]=vs[m][perm_norton]
  n_steps_norton[m]=n_steps_norton[m][perm_norton]
  error_bars_norton[m]=error_bars_norton[m][perm_norton]
  asymptotic_vars_norton[m]=asymptotic_vars_norton[m][perm_norton]
  
  scatter!(single_plot,dΛs[m],vs[m],markershape=:xcross,label="$(m)_N",markersize=2)
  scatter!(single_plot_linear_regime,dΛs[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_N",markersize=2)
  scatter!(joint_plot,dΛs[m],vs[m],markershape=:xcross,label="$(m)_N",markersize=2,xerror=error_bars_norton[m],msc=:auto)
  scatter!(joint_plot_linear_regime,dΛs[m][1:n_linear_regime],vs[m][1:n_linear_regime],markershape=:xcross,label="$(m)_N",markersize=2,xerror=error_bars_norton[m])
  scatter!(plot_asympt_var,dΛs[m],asymptotic_vars_norton[m],markershape=:xcross,markersize=2,label="$(m)_N")
  scatter!(plot_asympt_var_linear_regime,dΛs[m][1:n_linear_regime],asymptotic_vars_norton[m][1:n_linear_regime],markershape=:xcross,markersize=2,label="$(m)_N")
  scatter!(plot_nsteps,dΛs[m],n_steps_norton[m],markershape=:xcross,markersize=2,label="$(m)_N")
  savefig(single_plot,"$(m).pdf")
  savefig(single_plot_linear_regime,"$(m)_linear_regime.pdf")
end

savefig(joint_plot,"joint_plot.pdf")
savefig(joint_plot_linear_regime,"joint_plot_linear_regime.pdf")
savefig(plot_nsteps,"nsteps_plot.pdf")
savefig(plot_asympt_var,"asymptotic_vars.pdf")
savefig(plot_asympt_var_linear_regime,"asymptotic_vars_linear_regime.pdf")