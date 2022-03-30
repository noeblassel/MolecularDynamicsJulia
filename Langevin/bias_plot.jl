#!/usr/bin/env julia

using Plots,Statistics,CSV


Base.run(`./scp_bias_files.sh`)

data_folder="./bias_dumps/"

files=["BABO.csv","BAOAB.csv","BAOA.csv","BAO.csv"]#,"GHMC.csv"]

n_regr=3#number of regression points
n_start_regr=1
y_lim_tol=0.005
explosion_threshold=60.0
orders=Dict("BAO"=>1,"BAOA"=> 2,"BAOAB"=>2,"BABO"=>2,"GHMC"=>0)
colors=Dict("BAO"=>:blue,"BAOA"=>:orange,"BAOAB"=>:green,"BABO"=>:red,"GHMC"=>:purple)

X=zeros(n_regr*length(files),1+length(files))
V_Y=zeros(n_regr*length(files))
K_Y=zeros(n_regr*length(files))
W_Y=zeros(n_regr*length(files))
legend_pos=:bottomleft

plot_V=plot(xlabel="Δt",ylabel="Potential Energy")
plot_K=plot(xlabel="Δt",ylabel="Kinetic Energy")
plot_W=plot(xlabel="Δt",ylabel="Virial")

dts=nothing


for (j,input_file) in enumerate(files)
    scheme=split(input_file,'.')[1]
    dat=CSV.File(data_folder*input_file;header=false,types=Float64)
    dat=[r for r in dat if all(!isa(x,Missing) for x in r)&&(maximum(abs.(r))<explosion_threshold)]
    order= orders[scheme]    

    if j==1
        global dts=Set(r[1] for r in dat)
        global dts=[dt for dt in dts]
        sort!(dts)
        X[:,1].=1
    end
    for (i,Δt) in enumerate(dts[n_start_regr:n_start_regr+n_regr-1])
        X[(j-1)*n_regr+i,j+1]=Δt^order
    end

    V=[mean(r[2] for r in dat if r[1]==dt) for dt in dts]
    K=[mean(r[3] for r in dat if r[1]==dt) for dt in dts]
    W=[mean(r[4] for r in dat if r[1]==dt) for dt in dts]

    for i in 1:n_regr
        V_Y[(j-1)*n_regr+i]=V[i+n_start_regr-1]
        K_Y[(j-1)*n_regr+i]=K[i+n_start_regr-1]
        W_Y[(j-1)*n_regr+i]=W[i+n_start_regr-1]
    end
    
    color= colors[scheme]

    scatter!(plot_V,dts,V,label="",markershape=:xcross,color=color)
    scatter!(plot_K,dts,K,label="",markershape=:xcross,color=color)
    scatter!(plot_W,dts,W,label="",markershape=:xcross,color=color)

end


#solve least squares

H=inv(transpose(X)*X)*transpose(X)
θ_V=H*V_Y
θ_K=H*K_Y
θ_W=H*W_Y


for (j,input_file) in enumerate(files)
    scheme=split(input_file,'.')[1]
    order= orders[scheme]
    color= colors[scheme]

    f_V(x)=θ_V[1]+θ_V[1+j]*x^order
    f_K(x)=θ_K[1]+θ_K[1+j]*x^order
    f_W(x)=θ_W[1]+θ_W[1+j]*x^order

    plot!(plot_V,f_V,0,maximum(dts),label=scheme,color=color,linestyle=:dot,legend=legend_pos)
    plot!(plot_K,f_K,0,maximum(dts),label=scheme,color=color,linestyle=:dot,legend=legend_pos)
    plot!(plot_W,f_W,0,maximum(dts),label=scheme,color=color,linestyle=:dot,legend=legend_pos)

end

println("Estimated virial: $(θ_W[1])")
println("Estimated potential energy: $(θ_V[1])")
println("Estimated kinetic energy: $(θ_K[1])")

"""ylims!(plot_V,(θ_V[1]-y_lim_tol*abs(θ_V[1]),θ_V[1]+y_lim_tol*abs(θ_V[1])))
ylims!(plot_K,(θ_K[1]-y_lim_tol*abs(θ_K[1]),θ_K[1]+y_lim_tol*abs(θ_K[1])))
ylims!(plot_W,(θ_W[1]-y_lim_tol*abs(θ_W[1]),θ_W[1]+y_lim_tol*abs(θ_W[1])))"""

savefig(plot_V,"./potential_energy_bias.pdf")
savefig(plot_K,"./kinetic_energy_bias.pdf")
savefig(plot_W,"./virial_bias.pdf")
