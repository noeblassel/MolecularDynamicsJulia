#!/usr/bin/env julia
using Plots

run(`./scp_files.sh`)
methods=["SINUSOIDAL","CONSTANT","LINEAR"]
prefixes=["velocity_thevenin_","forcing_norton_"]
suffix=".out"

label_dict=Dict("velocity_thevenin_"=>"velocity","forcing_norton_"=>"forcing")
label_dict_opp=Dict("velocity_thevenin_"=>"forcing","forcing_norton_"=>"velocity")

for method in methods
    println(method)
    for prefix in prefixes
        if isfile(prefix*method*suffix)

            println("\t",prefix)
            f=open(prefix*method*suffix,"r")

            _,Ly=split(readline(f))
            Ly=parse(Float64,Ly)


            f_dict=Dict("SINUSOIDAL"=>(y::Float64 -> sin(2π*y/Ly)),"CONSTANT"=>(y::Float64 -> (y<Ly/2) ? -1 : 1),"LINEAR"=>(y::Float64 -> (y<Ly/2) ? 4*(y-Ly/4)/Ly : 4*(3Ly/4-y)/Ly))

            _,n_bins=split(readline(f))
            n_bins=parse(Int64,n_bins)

            _,n_samples=split(readline(f))
            n_samples=parse(Int64,n_samples)

            _,v=split(readline(f))
            v=parse(Float64,v)

            profile=parse.(Float64,split(readline(f)))/(v*n_samples)

            y_range=(0:n_bins-1)*(Ly/n_bins)

            F=f_dict[method].(y_range)

            plot(profile,y_range,xlabel="reduced magnitude",ylabel="y",label=label_dict[prefix],color=:blue,legend=:topright)
            plot!(F,y_range,label=label_dict_opp[prefix],linestyle=:dot,color=:red)

            savefig(prefix*method*".pdf")
        
        end
    end
end
            