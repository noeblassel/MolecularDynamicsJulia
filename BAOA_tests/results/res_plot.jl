#!/bin/env julia

include("../potential.jl")
using Plots
run(`./scp_files.sh`)

dts=[0.1,0.15,0.2,0.25,0.3,0.35,0.4]
methods=["BAOA","BAOAB"]
potentials=["PERIODIC","QUADRATIC","DOUBLE_WELL"]

h_margin=0.3
v_margin=0.2
w_side=1000.0

plotsize_joint=(w_side*(1+h_margin)*length(methods),w_side*(1+v_margin)*length(dts))
plotsize_marginal=(w_side*(1+h_margin)*2,w_side*(1+v_margin)*length(dts))

for potential in potentials
    println(potential)

    joint_plots=[]
    marginal_plots=[]

    p_lims=(-5.0,5.0)

    if potential=="PERIODIC"
        q_lims=(0.0,1.0)
    else
        q_lims=(-5.0,5.0)
    end

    κ(p)=exp(-0.5p^2)/sqrt(2π)
    ar=(last(q_lims)-first(q_lims))/(last(p_lims)-first(p_lims))
    L=1.0

    if potential=="PERIODIC"
        Z=1.2660648774739363
        V=periodic_potential
    elseif potential=="QUADRATIC"
        Z=sqrt(2π)
        V=quadratic_potential
    elseif potential=="DOUBLE_WELL"
        Z=1.0
        V=double_well_potential
    end

    ν(q)=exp(-V(q,L))/Z

    for dt in dts
        println("\t",dt)

        p_marginal_plot=plot(xlims=p_lims,xlabel="p",ylabel="κ",size=(w_side,w_side))
        q_marginal_plot=plot(x_lims=q_lims,xlabel="q",ylabel="ν",size=(w_side,w_side))

        plot!(p_marginal_plot,κ,label="reference")
        plot!(q_marginal_plot,ν,label="reference")

        joint_subplots=[]

        for method in methods
                println("\t\t",method)
                rows=split(read("bins_$(method)_$(potential)_$(dt).out",String),';')
                M=reduce(hcat,[parse.(Int64,split(r)) for r in rows[1:end-1]])
                sum_M=sum(M)
                (Np,Nq)=size(M)
                

                dp=(last(p_lims)-first(p_lims))/Np
                dq=(last(q_lims)-first(q_lims))/Nq

                D=M/(sum_M*dp*dq)

                prange=range(p_lims...,Np)
                qrange=range(q_lims...,Nq)

                push!(joint_subplots,heatmap(qrange,prange,D,c=:hsv,xlabel="q",ylabel="p",aspect_ratio=ar,title="$(method)",size=(w_side,w_side)))#,clims=(0.0,1.0)

                q_marginal=zeros(Nq)
                p_marginal=zeros(Np)

                for i=1:Nq
                    I=0
                    for (p1,p2)=zip(D[:,i],D[2:end,i])
                        I+=dp*(p1+p2)/2
                    end
                    q_marginal[i]=I
                end

                for i=1:Np
                    I=0
                    for (q1,q2)=zip(D[i,:],D[i,2:end])
                        I+=dq*(q1+q2)/2
                    end
                    p_marginal[i]=I
                end

                plot!(p_marginal_plot,prange,p_marginal,label=method,linestyle=:dot)
                plot!(q_marginal_plot,qrange,q_marginal,label=method,linestyle=:dot)
        end
        push!(joint_plots,plot(joint_subplots...,layout=(1,length(methods)),plot_title="Δt=$(dt)"))
        push!(marginal_plots,plot(p_marginal_plot,q_marginal_plot,layout=(1,2),plot_title="Δt=$(dt)"))
        savefig(last(joint_plots),"./plots/joint/joint_$(potential)_$(dt).pdf")
        savefig(last(marginal_plots),"./plots/marginal/marginal_$(potential)_$(dt).pdf")
    end
    plot(joint_plots...,layout=(length(dts),1),size=plotsize_joint)
    savefig("./plots/joint_distributions_$(potential).pdf")
    plot(marginal_plots...,layout=(length(dts),1),size=plotsize_marginal)
    savefig("./plots/marginal_distributions_$(potential).pdf")
end
