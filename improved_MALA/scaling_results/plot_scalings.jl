using Plots, LinearAlgebra, Base.Meta

N_start=12
N_regr=4

Dict(("hmc")=>:green,("em")=>:blue,()=>:red)
normal_plot=plot(xlabel="Δt",ylabel="rejection rate",legend=:bottomright,xlims=(0,1e-3),ylims=(0.0,0.6))
log_plot=plot(xlabel="Δt",ylabel="rejection rate",xaxis=:log,yaxis=:log,legend=:bottomright)

col=Dict("em"=>:red,"hmc"=>:blue)

for rule in ["metropolis","barker"]
    if rule=="barker"
        f=open("barker_hmc_scaling.out")
        s=readline(f)
        eval(Meta.parse(s))
        s=readline(f)
        eval(Meta.parse(s))
        s=readline(f)
        eval(Meta.parse(s))
        close(f)

        lg_dts=log.(dts[N_start:N_start+N_regr])
        lg_dts .-=first(lg_dts)

        lg_R_rel=log.(R_rel[N_start:N_start+N_regr])
        lg_R_rel .-= first(lg_R_rel)
        α_rel=inv(dot(lg_dts,lg_dts))*dot(lg_dts,lg_R_rel)
        f_regr_rel(t)= R_rel[N_start]*(t^α_rel)/(dts[N_start]^α_rel)

        #dts=dts[N_start:end]
        #R_abs=R_abs[N_start:end]
        #R_rel=R_rel[N_start:end]

        scatter!(normal_plot,dts[N_start:end],R_rel[N_start:end],color=:green,markershape=:xcross,label="barker hmc avg",xlims=(0,maximum(dts)))
        plot!(normal_plot,f_regr_rel,color=:green,linewidth=0.5,linestyle=:dash,label="")

        scatter!(log_plot,dts[N_start:end],R_rel[N_start:end],color=:green,markershape=:xcross,label="barker hmc avg")
        plot!(log_plot,f_regr_rel,color=:green,linewidth=0.5,linestyle=:dash,label="slope: $(round(α_rel,digits=2))")
        
        lg_R_abs=log.(R_abs[N_start:N_start+N_regr])
        lg_R_abs .-= first(lg_R_abs)
        α_abs=inv(dot(lg_dts,lg_dts))*dot(lg_dts,lg_R_abs)
        f_regr_abs(t)= R_abs[N_start]*(t^α_abs)/(dts[N_start]^α_abs)

        scatter!(normal_plot,dts[N_start:end],R_abs[N_start:end],color=:purple,markershape=:xcross,label="barker hmc abs",xlims=(0,maximum(dts)))
        plot!(normal_plot,f_regr_abs,color=:purple,linewidth=0.5,linestyle=:dash,label="")

        scatter!(log_plot,dts[N_start:end],R_abs[N_start:end],color=:purple,markershape=:xcross,label="barker hmc abs")
        plot!(log_plot,f_regr_abs,color=:purple,linewidth=0.5,linestyle=:dash,label="slope: $(round(α_abs,digits=2))")

    else
        for proposal in ["em","hmc"]
            f=open("$(rule)_$(proposal)_scaling.out")
            s=readline(f)
            eval(Meta.parse(s))
            s=readline(f)
            eval(Meta.parse(s))
            close(f)

            lg_dts=log.(dts[N_start:N_start+N_regr])
            lg_dts .-=first(lg_dts)

            lg_R=log.(R[N_start:N_start+N_regr])
            lg_R .-= first(lg_R)

            α=inv(dot(lg_dts,lg_dts))*dot(lg_dts,lg_R)
            f_regr(t)= R[N_start]*(t^α)/(dts[N_start]^α)

            #dts=dts[N_start:end]
            #R=R[N_start:end]

            scatter!(normal_plot,dts[N_start:end],R[N_start:end],color=col[proposal],markershape=:xcross,label="$(rule) $(proposal)",xlims=(0,maximum(dts)))
            plot!(normal_plot,f_regr,color=col[proposal],linewidth=0.5,linestyle=:dash,label="")

            scatter!(log_plot,dts[N_start:end],R[N_start:end],color=col[proposal],markershape=:xcross,label="$(rule) $(proposal)")
            plot!(log_plot,f_regr,color=col[proposal],linewidth=0.5,linestyle=:dash,label="slope: $(round(α,digits=2))")

        end

    end

end

savefig(normal_plot,"scalings.pdf")
savefig(log_plot,"scalings_log.pdf")