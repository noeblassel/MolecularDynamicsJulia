using Plots,LinearAlgebra

include("potential.jl")
dts=[0.02,0.04,0.06,0.08,0.1,0.12,0.14]
methods=["BAOA","BAOAB"]
potential=ARGS[1]

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential,"TILTED_DOUBLE_WELL"=>tilted_double_well_potential)
force_dict=Dict("PERIODIC"=>minus_d_periodic_potential,"QUADRATIC"=>minus_d_quadratic_potential,"DOUBLE_WELL"=>minus_d_double_well_potential,"TILTED_DOUBLE_WELL"=>minus_d_tilted_double_well_potential)

qlims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-5.0,5.0),"DOUBLE_WELL"=>(-5.0,5.0),"TILTED_DOUBLE_WELL"=>(-5.0,5.0))
plims=(-5.0,5.0)

L=1.0
F=force_dict[potential]
qlims=qlims_dict[potential]
phi(q,p)=-p*F(q,L)#error should be of order 1 on this observable 



normal_plot=plot(xlims=(minimum(dts),maximum(dts)),xlabel="Δt",ylabel="E[g]",legend=:outertopright)
ll_plot=plot(xaxis=:log,yaxis=:log,xlims=(minimum(dts),maximum(dts)),xlabel="Δt",ylabel="E[g]",legend=:outertopright)

for method in methods
    println(method)
    Is=[]
    for dt=dts
        println("\t",dt)
        rows=readlines("results/bins_$(method)_$(potential)_$(dt).out")
        M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
        sum_M=sum(M)
        (Nq,Np)=size(M)
        
        dp=(last(plims)-first(plims))/Np
        dq=(last(qlims)-first(qlims))/Nq

        D=M/(sum_M*dq*dp)

        prange=range(plims...,Np)
        qrange=range(qlims...,Nq)

        I_p=zero(qrange)

        for (i,q)=enumerate(qrange)
            I=0
            for j=eachindex(prange[1:end-1])
                I+=dp*(phi(q,prange[j])*D[i,j]+phi(q,prange[j+1])*D[i,j+1])/2
            end
            push!(I_p,I)
        end

        I=0
        for (I_p1,I_p2)=zip(I_p,I_p[2:end])
            I+=(I_p1+I_p2)*dq/2
        end

        push!(Is,I)
    end
    a=inv(dot(dts,dts))*dot(Is,dts)
    scatter!(ll_plot,dts,abs.(Is),markershape=:xcross,label="$(method)")
    scatter!(normal_plot,dts,abs.(Is),markershape=:xcross,label="$(method)")
    #plot!(t->a*t,0,maximum(dts),linestyle=:dot,color=:red,label="first order fit")

end

plot!(ll_plot,x->x^4,minimum(dts),maximum(dts),linestyle=:dot,label="Δt⁴")
plot!(ll_plot,x->x,minimum(dts),maximum(dts),linestyle=:dot,label="Δt")
savefig(ll_plot,"$(potential)_bias_loglog.pdf")
savefig(normal_plot,"$(potential)_bias.pdf")