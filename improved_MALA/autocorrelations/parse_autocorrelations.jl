using Plots

file_prefix="autocorrelation_"
file_bodies=["BARKER_HMC","METROPOLIS_HMC","METROPOLIS_EM","BARKER_EM"]
file_suffix=".out"

dt=5e-4
N=64
D=3

norm_plot=plot(xlabel="t",ylabel="∇V/√dN autocorrelation")
log_plot=plot(xlabel="t",ylabel="∇V/√dN autocorrelation",yaxis=:log)

N_plot=1000
for bod in file_bodies
    dump=open(file_prefix*bod*file_suffix,"r")
    n_vecs=read(dump,Int64)

    for i=1:6*n_vecs#discard average of gradient
        read(dump,Float64)
    end

    sum_sq_A=read(dump,Float64)
    read(dump,Float64)#same

    n_corr=read(dump,Int64)
    corr=Float64[]
    for i=1:n_corr
        push!(corr,read(dump,Float64))
    end

    n_samples=read(dump,Int64)
    println(n_samples)
    for i=1:n_corr
        corr[i]#/=(N*D)*(n_samples-i+1)
    end
    println(sum_sq_A," ",corr[1])
    dt_range=dt*(0:(n_corr-1))
    plot!(norm_plot,dt_range[1:N_plot],corr[1:N_plot],label=bod,linewidth=0.5)
    plot!(log_plot,dt_range[2:N_plot],abs.(corr[2:N_plot]),yaxis=:log,label=bod,linewidth=0.5)
    
end

savefig(log_plot,"autocorrelations_logy.pdf")
savefig(norm_plot,"autocorrelations.pdf")