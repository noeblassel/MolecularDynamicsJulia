using Plots, LinearAlgebra

cargs=["clustern15:/libre/blasseln/MolecularDynamicsJulia/improved_MALA/1d/*.out","."]
run(`scp $cargs`)
files=readdir()

n_start=11
n_plot=30

n_regr_start=36
n_regr=5

dts_me=Float64[]
gk_me=Float64[]
v_me=Float64[]
e_me=Float64[]

dts_mh=Float64[]
e_mh=Float64[]
gk_mh=Float64[]
v_mh=Float64[]

dts_bh=Float64[]
e_bh=Float64[]
gk_bh=Float64[]
v_bh=Float64[]

Diff=0.62386037

for f in files

    if endswith(f,".out") && (f!="nohup.out")
        (rule,proposal,l)=match(r"^(.+)_(.+)_(.+)\.out",f).captures
        lg_dt=parse(Float64,l)
        dt=10^lg_dt
        lines=readlines(f)
        N=length(lines)
        C0_sum=0.0
        I_sum=0.0
        D_sum=0.0
        for l in lines
            (C0,I,D)=parse.(Float64,split(l))
            C0_sum+=C0
            I_sum+=I
            D_sum+=D
        end
        C0hat=C0_sum/N
        Ihat=I_sum/N
        Dhat=D_sum/N

        if (rule,proposal)==("metropolis","EM")
            push!(gk_me,1-Ihat)
            push!(e_me,Dhat)
            push!(dts_me,dt)
            push!(v_me,C0hat)
        elseif (rule,proposal)==("metropolis","HMC")
            push!(dts_mh,dt)
            push!(e_mh,Dhat)
            push!(gk_mh,1-(Ihat-dt*C0hat/2))
            push!(v_mh,C0hat)
        else
            push!(dts_bh,dt)
            push!(e_bh,2*Dhat)
            push!(gk_bh,1-(Ihat-dt*C0hat)/2)
            push!(v_bh,C0hat)

        end
    end

end
normal_plot=plot(legend=:bottomleft,xlims=(0.0,1e-2),xlabel="Δt",ylabel="estimated diffusion coefficient")
log_plot=plot(legend=:bottomright,xaxis=:log,yaxis=:log,xlabel="Δt",ylabel="absolute error")
eq_plot=plot()

perm=sortperm(dts_me)
dts_me=dts_me[perm]
gk_me=gk_me[perm]
abs_gk_me=abs.(gk_me .- Diff)
e_me=e_me[perm]
abs_e_me=abs.(e_me .- Diff)
v_me=v_me[perm]

lg_dts_me=log.(dts_me[n_regr_start:n_regr_start+n_regr]).-log(dts_me[n_regr_start])
lg_gk_me=log.(abs_gk_me[n_regr_start:n_regr_start+n_regr]).-log(abs_gk_me[n_regr_start])
lg_e_me=log.(abs_e_me[n_regr_start:n_regr_start+n_regr]).-log(abs_e_me[n_regr_start])

a_gk_me=inv(dot(lg_dts_me,lg_dts_me))*dot(lg_dts_me,lg_gk_me)
a_e_me=inv(dot(lg_dts_me,lg_dts_me))*dot(lg_dts_me,lg_e_me)

scatter!(normal_plot,dts_me[n_start:n_start+n_plot],gk_me[n_start:n_start+n_plot],label="gk-me",markershape=:xcross,markersize=3)
scatter!(normal_plot,dts_me[n_start:n_start+n_plot],e_me[n_start:n_start+n_plot],label="e-me",markershape=:xcross,markersize=3)

scatter!(eq_plot,dts_me,v_me,markershape=:xcross,label="me")

scatter!(log_plot,dts_me[n_start:n_start+n_plot],abs_gk_me[n_start:n_start+n_plot],label="gk-me($(round(a_gk_me,digits=2)))",markershape=:xcross,markersize=3,color=:red)
scatter!(log_plot,dts_me[n_start:n_start+n_plot],abs_e_me[n_start:n_start+n_plot],label="e-me($(round(a_e_me,digits=2)))",markershape=:xcross,markersize=3,color=:blue)
plot!(log_plot,t->t^a_gk_me*(abs_gk_me[n_regr_start]/dts_me[n_regr_start]^a_gk_me),linestyle=:dot,label="",color=:red)
plot!(log_plot,t->t^a_e_me*(abs_e_me[n_regr_start]/dts_me[n_regr_start]^a_e_me),linestyle=:dot,label="",color=:blue)

perm=sortperm(dts_mh)
dts_mh=dts_mh[perm]
gk_mh=gk_mh[perm]
e_mh=e_mh[perm]
v_mh=v_mh[perm]
abs_gk_mh=abs.(gk_mh .- Diff)
abs_e_mh=abs.(e_mh .- Diff)

lg_dts_mh=log.(dts_mh[n_regr_start:n_regr_start+n_regr]).-log(dts_mh[n_regr_start])
lg_gk_mh=log.(abs_gk_mh[n_regr_start:n_regr_start+n_regr]).-log(abs_gk_mh[n_regr_start])
lg_e_mh=log.(abs_e_mh[n_regr_start:n_regr_start+n_regr]).-log(abs_e_mh[n_regr_start])

a_gk_mh=inv(dot(lg_dts_mh,lg_dts_mh))*dot(lg_dts_mh,lg_gk_mh)
a_e_mh=inv(dot(lg_dts_mh,lg_dts_mh))*dot(lg_dts_mh,lg_e_mh)

scatter!(normal_plot,dts_mh[n_start:n_start+n_plot],gk_mh[n_start:n_start+n_plot],label="gk-mh",markershape=:xcross,markersize=3,color=:green)
scatter!(normal_plot,dts_mh[n_start:n_start+n_plot],e_mh[n_start:n_start+n_plot],label="e-mh",markershape=:xcross,markersize=3,color=:purple)

scatter!(eq_plot,dts_mh,v_mh,markershape=:xcross,label="mh")

scatter!(log_plot,dts_mh[n_start:n_start+n_plot],abs_gk_mh[n_start:n_start+n_plot],label="gk-mh($(round(a_gk_mh,digits=2)))",markershape=:xcross,markersize=3,color=:green)
scatter!(log_plot,dts_mh[n_start:n_start+n_plot],abs_e_mh[n_start:n_start+n_plot],label="e-mh($(round(a_e_mh,digits=2)))",markershape=:xcross,markersize=3,color=:purple)
plot!(log_plot,t->t^a_gk_mh*(abs_gk_mh[n_regr_start]/dts_mh[n_regr_start]^a_gk_mh),linestyle=:dot,label="",color=:green)
plot!(log_plot,t->t^a_e_mh*(abs_e_mh[n_regr_start]/dts_mh[n_regr_start]^a_e_mh),linestyle=:dot,label="",color=:purple)

perm=sortperm(dts_bh)
dts_bh=dts_bh[perm]
gk_bh=gk_bh[perm]
e_bh=e_bh[perm]
v_bh=v_bh[perm]
abs_gk_bh=abs.(gk_bh .- Diff)
abs_e_bh=abs.(e_bh .- Diff)
scatter!(normal_plot,dts_bh[n_start:n_start+n_plot],gk_bh[n_start:n_start+n_plot],label="gk-bh",markershape=:xcross,markersize=3)
scatter!(normal_plot,dts_bh[n_start:n_start+n_plot],e_bh[n_start:n_start+n_plot],label="e-bh",markershape=:xcross,markersize=3)
#scatter!(log_plot,dts_bh[n_start:n_start+n_plot],abs_gk_bh[n_start:n_start+n_plot],label="gk-bh",markershape=:xcross,markersize=3)
#scatter!(log_plot,dts_bh[n_start:n_start+n_plot],abs_e_bh[n_start:n_start+n_plot],label="e-bh",markershape=:xcross,markersize=3)
scatter!(eq_plot,dts_bh,v_bh,markershape=:xcross,label="bh")
#plot!(x->x,linestyle=:dot,label="1")
#plot!(x->x^1.5,linestyle=:dot,label="1.5")
#plot!(x->x^2,linestyle=:dot,label="2")
savefig(normal_plot,"plots/D.pdf")
savefig(log_plot,"plots/D_log.pdf")
savefig(eq_plot,"plots/eq.pdf")