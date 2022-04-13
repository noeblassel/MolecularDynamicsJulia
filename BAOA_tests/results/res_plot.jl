#!/bin/env julia

include("../potential.jl")
using Plots
run(`./scp_files.sh`)

dts=[0.1,0.15,0.2,0.25,0.3,0.35,0.4]
methods=["BAOA","BAOAB"]
potentials=["PERIODIC","QUADRATIC","DOUBLE_WELL"]

N_ref_pts=1000

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential,"TILTED_DOUBLE_WELL"=>tilted_double_well_potential)
force_dict=Dict("PERIODIC"=>minus_d_periodic_potential,"QUADRATIC"=>minus_d_quadratic_potential,"DOUBLE_WELL"=>minus_d_double_well_potential,"TILTED_DOUBLE_WELL"=>minus_d_tilted_double_well_potential)
Z_dict=Dict("PERIODIC"=>1.2660648774739363,"QUADRATIC"=>sqrt(2π),"DOUBLE_WELL"=>0.9423576216810029,"TILTED_DOUBLE_WELL"=>2.3080545166317865) #numerically computed using a trapezoid rule with dq=1e-6

qlims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-5.0,5.0),"DOUBLE_WELL"=>(-5.0,5.0),"TILTED_DOUBLE_WELL"=>(-5.0,5.0))
plims=(-5.0,5.0)

κ(p)=exp(-0.5p^2)/sqrt(2π)
L=1.0


color_lims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(0.0,0.2),"DOUBLE_WELL"=>(0.0,0.25),"TILTED_DOUBLE_WELL"=>(0.0,0.25))

h_margin=0.3
v_margin=0.2
w_side=1000.0
aspect_ratio=1.41

plotsize_joint=(w_side*(1+h_margin)*(length(methods)+1),w_side*(1+v_margin)*length(dts)/aspect_ratio)
plotsize_marginal=(w_side*(1+h_margin)*2,w_side*(1+v_margin)*length(dts)/aspect_ratio)

for potential in potentials
    
    println(potential)
    V=potential_dict[potential]
    F=force_dict[potential]
    ν(q)=exp(-V(q,L))/Z_dict[potential]
    μ(q,p)=ν(q)*κ(p)
    first_order_term(q,p)=-p*F(q,L)*μ(q,p)/2#from TU lemma knowing BAOAB is second order + Taylor expansion 
    qlims=qlims_dict[potential]

    MU=zeros(N_ref_pts,N_ref_pts)
    FOT=zeros(N_ref_pts,N_ref_pts)

    ref_prange=range(plims...,N_ref_pts)
    ref_qrange=range(qlims...,N_ref_pts)
    
    #compute reference terms

    for i=eachindex(MU)
        iq=1+mod(i-1,N_ref_pts)
        ip=1+div(i-1,N_ref_pts)
        q=ref_qrange[iq]
        p=ref_prange[ip]
        MU[i]=μ(q,p)
        FOT[i]=first_order_term(q,p)
    end
    
    marginal_plots=[]
    joint_plots=[]

    ar=(last(qlims)-first(qlims))/(last(plims)-first(plims))    

    for dt in dts
        println("\t",dt)

        p_marginal_plot=plot(xlims=plims,xlabel="p",ylabel="κ",size=(w_side,w_side/aspect_ratio),legend=:topleft)
        q_marginal_plot=plot(x_lims=qlims,xlabel="q",ylabel="ν",size=(w_side,w_side/aspect_ratio),legend=:topleft)

        plot!(p_marginal_plot,κ,label="reference")
        plot!(q_marginal_plot,ν,label="reference")

        joint_subplots=[]

        for method in methods
                println("\t\t",method)
                rows=readlines("bins_$(method)_$(potential)_$(dt).out")
                M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
                sum_M=sum(M)
                (Nq,Np)=size(M)
                

                dp=(last(plims)-first(plims))/Np
                dq=(last(qlims)-first(qlims))/Nq

                D=M/(sum_M*dp*dq)

                prange=range(plims...,Np)
                qrange=range(qlims...,Nq)

                push!(joint_subplots,heatmap(qrange,prange,D,c=:hsv,xlabel="q",ylabel="p",aspect_ratio=ar,title="$(method)",size=(w_side,w_side/aspect_ratio),clims=color_lims_dict[potential]))#

                q_marginal=zeros(Nq)
                p_marginal=zeros(Np)

                #compute approximate marginal distributions --- trapezoid rule ---

                for i=1:Nq
                    I=0
                    for (p1,p2)=zip(D[i,:],D[i,2:end])
                        I+=dp*(p1+p2)/2
                    end
                    q_marginal[i]=I
                end

                for i=1:Np
                    I=0
                    for (q1,q2)=zip(D[:,i],D[2:end,i])
                        I+=dq*(q1+q2)/2
                    end
                    p_marginal[i]=I
                end

                plot!(p_marginal_plot,prange,p_marginal,label=method,linestyle=:dot)
                plot!(q_marginal_plot,qrange,q_marginal,label=method,linestyle=:dot)
        end
        push!(joint_subplots,heatmap(ref_qrange,ref_prange,MU+dt*FOT,c=:hsv,xlabel="q",ylabel="p",aspect_ratio=ar,title="theoretical BAOA first order",size=(w_side,w_side/aspect_ratio),clims=color_lims_dict[potential]))

        push!(joint_plots,plot(joint_subplots...,layout=(1,length(methods)+1),plot_title="Δt=$(dt)"))
        push!(marginal_plots,plot(p_marginal_plot,q_marginal_plot,layout=(1,2),plot_title="Δt=$(dt)"))
        savefig(last(joint_plots),"./plots/joint/joint_$(potential)_$(dt).pdf")
        savefig(last(marginal_plots),"./plots/marginal/marginal_$(potential)_$(dt).pdf")
    end
    plot(joint_plots...,layout=(length(dts),1),size=plotsize_joint)
    savefig("./plots/joint_distributions_$(potential).pdf")
    plot(marginal_plots...,layout=(length(dts),1),size=plotsize_marginal)
    savefig("./plots/marginal_distributions_$(potential).pdf")
end
