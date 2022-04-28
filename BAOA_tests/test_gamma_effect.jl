#!/bin/env julia

include("./potential.jl")
using Plots
#run(`./scp_files.sh`)


potential="DOUBLE_WELL"
dt=0.3

gammas=[0.1,1.0,10.0]
methods=["BAOAB","BAOA"]

N_ref_pts=1000

potential_dict=Dict("PERIODIC"=>periodic_potential,"QUADRATIC"=>quadratic_potential,"DOUBLE_WELL"=>double_well_potential,"TILTED_DOUBLE_WELL"=>tilted_double_well_potential)
plims=(-5.0,5.0)
qlims_dict=Dict("PERIODIC"=>(0.0,1.0),"QUADRATIC"=>(-5.0,5.0),"DOUBLE_WELL"=>(-5.0,5.0),"TILTED_DOUBLE_WELL"=>(-5.0,5.0))
qlims=qlims_dict[potential]

κ(p)=exp(-0.5p^2)/sqrt(2π)
L=1.0

linestyle_dict=Dict("BAOA"=>:dash,"BAOAB"=>:dot)


h_margin=0.3
v_margin=0.2
w_side=1000.0
aspect_ratio=1.41

plot(xlims=plims,xlabel="p",ylabel="κ",legend=:topleft)
plot!(κ,label="reference")

for γ in gammas
    
    println(γ)    

    ar=(last(plims)-first(plims))/(last(qlims)-first(qlims))

        

        joint_subplots=[]

        for method in methods
                println("\t\t",method)
                rows=readlines("./results/bins_$(method)_$(potential)_$(dt)_$(γ).out")
                M=reduce(vcat,[transpose(parse.(Int64,split(r))) for r in rows])
                sum_M=sum(M)
                (Nq,Np)=size(M)
                

                dp=(last(plims)-first(plims))/Np
                dq=(last(qlims)-first(qlims))/Nq

                D=M/(sum_M*dp*dq)

                prange=range(plims...,Np)
                qrange=range(qlims...,Nq)


                p_marginal=zeros(Np)

                #compute approximate marginal distributions --- trapezoid rule ---

                for i=1:Np
                    I=0
                    for (q1,q2)=zip(D[:,i],D[2:end,i])
                        I+=dq*(q1+q2)/2
                    end
                    p_marginal[i]=I
                end

                plot!(prange,p_marginal,label="$(method) γ=$(γ)",linestyle=linestyle_dict[method])
        end
        
end

savefig("gamma_effect_$(potential)_$(dt).pdf")