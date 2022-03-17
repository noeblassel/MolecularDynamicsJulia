#!/usr/bin/env julia

using Plots,Statistics,Polynomials


Base.run(`./scp_bias_files.sh`)

babo_file="babo.out"
baoab_file="baoab.out"
baoa_file="baoa.out"
bao_file="bao.out"

orders=Dict("bao"=>1,"baoa"=> 1,"baoab"=>2,"babo"=>2)
colors=Dict("bao"=>:blue,"baoa"=>:orange,"baoab"=>:green,"babo"=>:red)


plot_V=plot(xlabel="Δt",ylabel="Potential Energy")
plot_K=plot(xlabel="Δt",ylabel="Kinetic Energy")
plot_W=plot(xlabel="Δt",ylabel="Virial")

for input_file in [bao_file,babo_file,baoab_file,baoa_file]
    scheme=split(input_file,'.')[1]
    
    order= orders[scheme]
    s=read(input_file,String)
    A=split(s,"\n")
    A=[parse.((Float64,),split(l)) for l in A[3:end-1] if length(l)>10]
    A=reduce(hcat,A)'
    
    dts=Set(r[1] for r in eachrow(A))
    dts=[dt for dt in dts]
    V=[mean(r[2] for r in eachrow(A) if r[1]==dt) for dt in dts]
    K=[mean(r[3] for r in eachrow(A) if r[1]==dt) for dt in dts]
    W=[mean(r[4] for r in eachrow(A) if r[1]==dt) for dt in dts]
    
    #fits of the form α+βx^order
    X= dts .^ order
    
    P_V=fit(X,V,1)
    P_K=fit(X,K,1)
    P_W=fit(X,W,1)

    f_V(x)=P_V(x^order)
    f_K(x)=P_K(x^order)
    f_W(x)=P_W(x^order)
    
    color= colors[scheme]

    scatter!(plot_V,dts,V,label="",markershape=:xcross,color=color)
    plot!(plot_V,f_V,0,maximum(dts),label=scheme,linestyle=:dot,color=color,legend=:top)

    scatter!(plot_K,dts,K,label="",markershape=:xcross,color=color)
    plot!(plot_K,f_K,0,maximum(dts),label=scheme,linestyle=:dot,color=color,legend=:top)
    
    
    scatter!(plot_W,dts,W,label="",markershape=:xcross,color=color)
    plot!(plot_W,f_W,0,maximum(dts),label=scheme,linestyle=:dot,color=color,legend=:top)

    V_est=f_V(0)
    K_est=f_K(0)
    W_est=f_W(0)

    lg_dt=log.(dts)

    lg_V=log.(abs.(V.-V_est))
    lg_K=log.(abs.(K.-K_est))
    lg_W=log.(abs.(W.-W_est))

    pow_V=fit(lg_dt,lg_V,1)
    pow_K=fit(lg_dt,lg_K,1)
    pow_W=fit(lg_dt,lg_W,1)

    exp_V=pow_V(1)-pow_V(0)
    exp_K=pow_K(1)-pow_K(0)
    exp_W=pow_W(1)-pow_W(0)

    println("$(scheme) estimated orders: $(exp_V) (potential energy) $(exp_K) (kinetic energy) $(exp_W) (virial)")
    
end

savefig(plot_V,"potential_energy_bias.pdf")
savefig(plot_K,"kinetic_energy_bias.pdf")
savefig(plot_W,"virial_bias.pdf")

Base.run(`rm babo.out baoab.out baoa.out bao.out`)