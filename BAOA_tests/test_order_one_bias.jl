using Plots,LinearAlgebra

include("potential.jl")
dts=[0.04,0.06,0.08,0.1,0.12,0.14]
methods=["BAOA","BAOAB"]
potential=ARGS[1]

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential,"TILTED_DOUBLE_WELL"=>tilted_double_well_potential)
force_dict=Dict("PERIODIC"=>minus_d_periodic_potential,"QUADRATIC"=>minus_d_quadratic_potential,"DOUBLE_WELL"=>minus_d_double_well_potential,"TILTED_DOUBLE_WELL"=>minus_d_tilted_double_well_potential)
order_dict=Dict("BAOA"=>1,"BAOAB"=>2,"OBABO"=>2)
qlims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-5.0,5.0),"DOUBLE_WELL"=>(-5.0,5.0),"TILTED_DOUBLE_WELL"=>(-5.0,5.0))
plims=(-5.0,5.0)

L=1.0
F=force_dict[potential]
qlims=qlims_dict[potential]
phi(q,p)=-p*F(q,L)#error should be of order 1 on this observable 

normal_plot_g=plot(xlims=(0,maximum(dts)+0.02),xlabel="Δt",ylabel="E[g]",legend=:topleft,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=16)
ll_plot_g=plot(xaxis=:log,yaxis=:log,xlims=(minimum(dts),maximum(dts)+0.02),xlabel="Δt",ylabel="E[g]",legend=:left,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=16,xticks=(dts,map(x->"$(x)",dts)))

normal_plot_p=plot(xlims=(0,maximum(dts)),xlabel="Δt",ylabel="E[p]",legend=:topleft,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=18)
ll_plot_p=plot(xaxis=:log,yaxis=:log,xlims=(minimum(dts),maximum(dts)),xlabel="Δt",ylabel="E[p]",legend=:topleft,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=16)


normal_plot_p2=plot(xlims=(0,maximum(dts)),xlabel="Δt",ylabel="E[p²]",legend=:topleft,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=16)
ll_plot_p2=plot(xaxis=:log,yaxis=:log,xlims=(minimum(dts),maximum(dts)+0.2),xlabel="Δt",ylabel="E[p²]",legend=:topleft,titlefontsize=18,tickfontsize=12,legendfontsize=14,labelfontsize=16)

for method in methods
    println(method)
    I_gs=Float64[]
    I_p2s=Float64[]
    I_ps=Float64[]

    for dt=dts
        println("\t",dt)
        file= (method=="OBABO") ? "results/bins_$(method)_$(potential)_$(dt)_1.0.out" : "results/bins_$(method)_$(potential)_$(dt).out" 
        rows=readlines(file)
        M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
        sum_M=sum(M)
        (Nq,Np)=size(M)
        
        dp=(last(plims)-first(plims))/Np
        dq=(last(qlims)-first(qlims))/Nq

        prange=range(plims...,Np)
        qrange=range(qlims...,Nq)

        I_p_g=zero(qrange)
        I_p_p2=zero(qrange)
        I_p_p=zero(qrange)

        for (i,q)=enumerate(qrange)
            I_g=0.0
            I_p2=0.0
            I_p=0.0

            for j=eachindex(prange[1:end-1])
                I_g+=(phi(q,prange[j])*M[i,j]+phi(q,prange[j+1])*M[i,j+1])/2
                I_p2+=((prange[j]^2-1)*M[i,j]+(prange[j+1]^2-1)*M[i,j+1])/2
                I_p+=(prange[j]*M[i,j]+prange[j+1]*M[i,j+1])/2
            end

            I_p_g[i]=I_g
            I_p_p2[i]=I_p2
            I_p_p[i]=I_p
        end

        I_g=0.0
        I_p2=0.0
        I_p=0.0

        for (I_1,I_2)=zip(I_p_g,I_p_g[2:end])
            I_g+=(I_1+I_2)/2
        end

        for (I_1,I_2)=zip(I_p_p2,I_p_p2[2:end])
            I_p2+=(I_1+I_2)/2
        end

        for (I_1,I_2)=zip(I_p_p,I_p_p[2:end])
            I_p+=(I_1+I_2)/2
        end



        push!(I_gs,I_g/sum_M)
        push!(I_p2s,I_p2/sum_M)
        push!(I_ps,I_p/sum_M)
    end

    order=order_dict[method]
    a=inv(dot(dts .^ order ,dts .^ order))*dot(I_gs,dts .^ order)

    regr=dt->a*dt^order
    scatter!(ll_plot_g,dts,abs.(I_gs),markershape=:xcross,label="$(method)")
    scatter!(normal_plot_g,dts,abs.(I_gs),markershape=:xcross,label="$(method)")
    plot!(normal_plot_g,regr,0,maximum(dts),linestyle=:dot,label="order $(order) fit")

    scatter!(ll_plot_p,dts,abs.(I_ps),markershape=:xcross,label="$(method)")
    scatter!(normal_plot_p,dts,abs.(I_ps),markershape=:xcross,label="$(method)")

    scatter!(ll_plot_p2,dts,abs.(I_p2s),markershape=:xcross,label="$(method)")
    scatter!(normal_plot_p2,dts,abs.(I_p2s),markershape=:xcross,label="$(method)")

    lg_Igs=log.(abs.(I_gs))
    lg_dts=log.(dts)

    lg_Igs=lg_Igs .- first(lg_Igs)
    lg_dts=lg_dts .- first(lg_dts)

    lg_a=inv(dot(lg_dts,lg_dts))*dot(lg_Igs,lg_dts)
    lg_regr=dt->(dt^lg_a)*first(abs.(I_gs))/(first(dts)^lg_a)
    plot!(ll_plot_g,lg_regr,minimum(dts),maximum(dts),linestyle=:dot,label="slope=$(round(lg_a,digits=2))")

end

savefig(ll_plot_g,"$(potential)_bias_loglog_g.pdf")
savefig(normal_plot_g,"$(potential)_bias_g.pdf")

savefig(ll_plot_p,"$(potential)_bias_loglog_p.pdf")
savefig(normal_plot_p,"$(potential)_bias_p.pdf")

savefig(ll_plot_p2,"$(potential)_bias_loglog_p2.pdf")
savefig(normal_plot_p2,"$(potential)_bias_p2.pdf")