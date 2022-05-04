#!/bin/env julia

include("../potential.jl")
using Plots, LinearAlgebra

dts=[0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18]#,0.08]
methods=["BAOA","BAOAB"]
potentials=["PERIODIC"]
plims=(-5.0,5.0)
N=2000
γ=1.0
m=2.0
m_inv=inv(m)

κ(qx,qy)=exp(-0.5(m_inv*qx^2+qy^2))/(2π*sqrt(m))

prange=range(plims...,N)
dp=(last(plims)-first(plims))/N
D_ref=zeros(N,N)

for i=1:N
    for j=1:N
        qx=first(plims)+0.5dp+(i-1)*dp
        qy=first(plims)+0.5dp+(j-1)*dp
        D_ref[i,j]=κ(qx,qy)
    end
end

#savefig(heatmap(prange,prange,D_ref,xlabel="q1",ylabel="q2",label="",c=:hsv,aspect_ratio=1.0,clims=(0.0,0.14)),"plots/kappa_reference.pdf")

for potential in potentials
    println("$(potential)")
    tv_plot=plot(xlims=(0,maximum(dts)),xlabel="Δt",ylabel="total variation",ylims=(0,0.2),legend=:topleft)
    tv_plot_log=plot(xlabel="Δt",ylabel="total variation",xaxis=:log,yaxis=:log,xlims=(minimum(dts),maximum(dts)),legend=:left)
    
    for method in methods 
        println("\t$(method)")
        TVs=[]
        for dt in dts
            println("\t\t$(dt)")
            rows=readlines("bins_$(method)_$(potential)_$(dt)_$(γ).out")
            M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
            D=M/(sum(M)*dp^2)
            
            savefig(heatmap(prange,prange,D-D_ref,xlabel="q1",ylabel="q2",label="",c=:hsv,aspect_ratio=1.0,clims=(-0.014,0.014)),"plots/diff_$(method)_$(potential)_$(dt)_$(γ).pdf")
            savefig(heatmap(prange,prange,D,xlabel="q1",ylabel="q2",label="",c=:hsv,aspect_ratio=1.0,clims=(0,0.14)),"plots/kappa_$(method)_$(potential)_$(dt)_$(γ).pdf")
            tv=sum(abs.(D-D_ref))*dp^2
            push!(TVs,tv)
            
            f= open("tvs_$(method).txt","a")
            println(f,dt," ",tv)
            close(f)

        end
        scatter!(tv_plot,dts,TVs,label="$(method)",markershape=:xcross)
        scatter!(tv_plot_log,dts,TVs,label="$(method)",markershape=:xcross)

        log_tvs=log.(TVs)
        log_dts=log.(dts)

        log_tvs=log_tvs .- first(log_tvs)
        log_dts=log_dts .- first(log_dts)

        lg_a=inv(dot(log_dts,log_dts))*dot(log_tvs,log_dts)
        lg_regr=dt->(dt^lg_a)*first(TVs)/(first(dts)^lg_a)

        plot!(tv_plot_log,lg_regr,minimum(dts),maximum(dts),linestyle=:dot,label="slope=$(round(lg_a,digits=2))")
        plot!(tv_plot,lg_regr,0,maximum(dts),label="exponent=$(round(lg_a,digits=2))")
    end

    savefig(tv_plot,"plots/TVs_$(potential).pdf")
    savefig(tv_plot_log,"plots/TVs_log_$(potential).pdf")
end