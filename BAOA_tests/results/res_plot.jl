#!/bin/env julia

include("../potential.jl")
using Plots
run(`./scp_files.sh`)

dts=[0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]
methods=["BAOA","BAOAB"]
potentials=["PERIODIC"]#"QUADRATIC","DOUBLE_WELL","TILTED_DOUBLE_WELL"]

N_ref_pts=1000

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential,"TILTED_DOUBLE_WELL"=>tilted_double_well_potential)
force_dict=Dict("PERIODIC"=>grad_periodic_potential,"QUADRATIC"=>grad_quadratic_potential,"DOUBLE_WELL"=>grad_double_well_potential,"TILTED_DOUBLE_WELL"=>grad_tilted_double_well_potential)
d2_dict=Dict("PERIODIC"=>d2_periodic_potential,"QUADRATIC"=>d2_quadratic_potential,"DOUBLE_WELL"=>d2_double_well_potential,"TILTED_DOUBLE_WELL"=>d2_tilted_double_well_potential)
Z_dict=Dict("PERIODIC"=>1.2660648774739363,"QUADRATIC"=>sqrt(2π),"DOUBLE_WELL"=>0.9423576216810029,"TILTED_DOUBLE_WELL"=>2.3080545166317865) #numerically computed using a trapezoid rule with dq=1e-6

qlims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-5.0,5.0),"DOUBLE_WELL"=>(-5.0,5.0),"TILTED_DOUBLE_WELL"=>(-5.0,5.0))
plims=(-5.0,5.0)

q_lims_plot_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-3.0,3.0),"DOUBLE_WELL"=>(-3.0,3.0),"TILTED_DOUBLE_WELL"=>(-3.0,3.0))
p_lims_plot_dict=Dict("PERIODIC"=>(-3.0,3.0),"QUADRATIC"=>(-3.0,3.0),"DOUBLE_WELL"=>(-3.0,3.0),"TILTED_DOUBLE_WELL"=>(-3.0,3.0))


κ(p)=exp(-abs(p^3)/3)/2.57579863371
L=1.0


color_lims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(0.0,0.2),"DOUBLE_WELL"=>(0.0,0.25),"TILTED_DOUBLE_WELL"=>(0.0,0.25))
linestyle_dict=Dict("BAOA"=>:dash,"BAOAB"=>:dot)


h_margin=0.3
v_margin=0.2
w_side=1000.0
aspect_ratio=1.41

plotsize_joint=(w_side*(1+h_margin)*(length(methods)+1),w_side*(1+v_margin)*length(dts)/aspect_ratio)
plotsize_marginal=(w_side*(1+h_margin)*2,w_side*(1+v_margin)*length(dts)/aspect_ratio)

for potential in potentials
    
    println(potential)
    V=potential_dict[potential]
    F=(q,L)->force_dict[potential](q,L)
    d2=d2_dict[potential]
    ν(q)=exp(V(q,L))/Z_dict[potential]
    μ(q,p)=ν(q)*κ(p)
    first_order_term(q,p)=-p^2*F(q,L)*μ(q,p)/2#from TU lemma knowing BAOAB is second order + Taylor expansion 
    qlims=qlims_dict[potential]

    MU=zeros(N_ref_pts,N_ref_pts)
    FOT=zeros(N_ref_pts,N_ref_pts)

    ref_prange=range(plims...,N_ref_pts)
    ref_qrange=range(qlims...,N_ref_pts)
    
    qlims_plot=q_lims_plot_dict[potential]
    plims_plot=p_lims_plot_dict[potential]

    (qm_pl,qM_pl)=qlims_plot
    (pm_pl,pM_pl)=plims_plot

    ix_qm_pl=max(1,floor(Int32,N_ref_pts*(qm_pl-first(qlims))/(last(qlims)-first(qlims))))
    ix_qM_pl=min(N_ref_pts,floor(Int32,N_ref_pts*(qM_pl-first(qlims))/(last(qlims)-first(qlims))))
    ix_pm_pl=max(1,floor(Int32,N_ref_pts*(pm_pl-first(plims))/(last(plims)-first(plims))))
    ix_pM_pl=min(N_ref_pts,floor(Int32,N_ref_pts*(pM_pl-first(plims))/(last(plims)-first(plims))))

    prange_plot=range(plims_plot...,ix_pM_pl-ix_pm_pl+1)
    qrange_plot=range(qlims_plot...,ix_qM_pl-ix_qm_pl+1)

    #compute reference terms

    for i=eachindex(MU)
        iq=1+mod(i-1,N_ref_pts)
        ip=1+div(i-1,N_ref_pts)
        q=ref_qrange[iq]
        p=ref_prange[ip]
        MU[i]=μ(q,p)
        FOT[i]=first_order_term(q,p)
    end

    correction_term_kinetic_BAOAB=zeros(N_ref_pts)
    
    ref_dq=(last(ref_qrange)-first(ref_qrange))/N_ref_pts

    for i=1:N_ref_pts
        p=ref_prange[i]
        I=0
        for (q1,q2)=zip(ref_qrange,ref_qrange[2:end])
            I+=((p^2*d2(q1,L)-F(q1,L)^2)*μ(q1,p)+(p^2*d2(q2,L)-F(q2,L)^2)*μ(q2,p))*ref_dq/2
        end
        correction_term_kinetic_BAOAB[i]=I
    end

    plot(ref_prange,correction_term_kinetic_BAOAB)
    savefig("correction_term_kinetic_BAOAB_$(potential).pdf")


    
    marginal_plots=[]
    joint_plots=[]

    ar=(last(plims_plot)-first(plims_plot))/(last(qlims_plot)-first(qlims_plot))

    for dt in dts
        println("\t",dt)

        p_marginal_plot=plot(xlims=plims,xlabel="p",ylabel="κ",legend=:topleft,legendfontsize=14,labelfontsize=14)
        q_marginal_plot=plot(x_lims=qlims,xlabel="q",ylabel="ν",legend=:topleft,legendfontsize=14,labelfontsize=14)

        plot!(p_marginal_plot,κ,label="reference")
        plot!(q_marginal_plot,ν,label="reference")

        joint_subplots=[]

        for method in methods
                println("\t\t",method)
                rows=readlines("bins_$(method)_$(potential)_$(dt)_1.0.out")
                M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
                sum_M=sum(M)
                (Nq,Np)=size(M)

                dp=(last(plims)-first(plims))/Np
                dq=(last(qlims)-first(qlims))/Nq
                
                D=M/(sum_M*dp*dq)

                prange=range(plims...,Np)
                qrange=range(qlims...,Nq)

                savefig(heatmap(prange_plot,qrange_plot,D[ix_qm_pl:ix_qM_pl,ix_pm_pl:ix_pM_pl],c=:hsv,xlabel="p",ylabel="q",aspect_ratio=ar,title="$(method)",clims=color_lims_dict[potential]),"./plots/joint/joint_$(method)_$(potential)_$(dt).pdf")#

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

                plot!(p_marginal_plot,prange,p_marginal,label=method,linestyle=linestyle_dict[method])

                if method=="BAOAB"
                    plot!(p_marginal_plot,prange,p_marginal+dt^2*correction_term_kinetic_BAOAB/8,label="correction",linestyle=linestyle_dict["BAOAB"])
                end
                plot!(q_marginal_plot,qrange,q_marginal,label=method,linestyle=linestyle_dict[method])
        end
        savefig(heatmap(prange_plot,qrange_plot,(MU+dt*FOT)[ix_qm_pl:ix_qM_pl,ix_pm_pl:ix_pM_pl],c=:hsv,xlabel="p",ylabel="q",aspect_ratio=ar,clims=color_lims_dict[potential],title="BAOA first order expansion"),"./plots/joint/joint_theoretical_$(potential)_$(dt).pdf")
        savefig(p_marginal_plot,"./plots/marginal/marginal_p_$(potential)_$(dt).pdf")
        savefig(q_marginal_plot,"./plots/marginal/marginal_q_$(potential)_$(dt).pdf")
"""        push!(joint_plots,plot(joint_subplots...,layout=(1,length(methods)+1),plot_title="Δt=$(dt)"))
        push!(marginal_plots,plot(p_marginal_plot,q_marginal_plot,layout=(1,2),plot_title="Δt=$(dt)"))
        savefig(last(joint_plots),"./plots/joint/joint_$(potential)_$(dt).pdf")
        savefig(last(marginal_plots),"./plots/marginal/marginal_$(potential)_$(dt).pdf")"""
    end
    savefig(heatmap(prange_plot,qrange_plot,MU[ix_qm_pl:ix_qM_pl,ix_pm_pl:ix_pM_pl],c=:hsv,xlabel="p",ylabel="q",aspect_ratio=ar,clims=color_lims_dict[potential],title="μ"),"./plots/joint/joint_reference_$(potential).pdf")
    #plot(joint_plots...,layout=(length(dts),1),size=plotsize_joint,plot_title="")
    #savefig("./plots/joint_distributions_$(potential).pdf")
   # plot(marginal_plots...,layout=(length(dts),1),size=plotsize_marginal)
    #savefig("./plots/marginal_distributions_$(potential).pdf")
end
