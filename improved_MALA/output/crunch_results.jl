using Plots

#run(`./scp_results.sh`)
rules = ["BARKER", "METROPOLIS"]
proposals = ["HMC", "EM"]

file_prefix = "autocorrelation_"
file_bodies = ["BARKER_HMC", "METROPOLIS_HMC", "METROPOLIS_EM", "BARKER_EM"]
file_suffix = ".out"

"""
einstein=Dict((rule,proposal)=>Float64[] for rule=rules for proposal=proposals)
gk=Dict((rule,proposal)=>Float64[] for rule=rules for proposal=proposals)"""

N = 64
β = 1.0
d = 3
T_corr = 0.4


lg_dts = [-5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0]

dts= 10 .^ lg_dts
normal_plot=plot(xlims=(0,maximum(dts)),xlabel="Δt",ylabel="diffusion coefficient",legend=:outerleft,legendfontsize=6,size=(1000,600))
log_plot=plot(xaxis=:log,yaxis=:log,xlabel="Δt",ylabel="diffusion coefficient",legend=:outerleft,legendfontsize=6,size=(1000,600))
eq_plot=plot(xlims=(0,maximum(dts)),xlabel="Δt",ylabel="|∇V|²",legend=:outerleft,legendfontsize=6,size=(1000,600))

for rule in rules
    println(rule)
    for proposal in proposals

        ((rule,proposal)==("BARKER","EM")) && continue

        println("\t", proposal)
        dts_gk = Float64[]
        D_gk = Float64[]
        dts_einstein = Float64[]
        D_einstein = Float64[]
        grad_V=Float64[]
        for lg_dt in lg_dts
            println("\t\t", lg_dt)

            sd_file = false
            gk_file = false

            if isfile("autocorrelation_$(rule)_$(proposal)_$(lg_dt).out")
                gk_file = true
                push!(dts_gk, 10^lg_dt)
                dt = 10^lg_dt
                f = open("autocorrelation_$(rule)_$(proposal)_$(lg_dt).out", "r")

                samp_count = 0
                est_Ds = Float64[]
                est_grad_V=Float64[]
                while (!eof(f))
                    samp_count += 1
                    N_grads = read(f, Int64)

                    #discard average gradient
                    for i = 1:2*d*N_grads
                        read(f, Float64)
                    end

                    sum_sq_grad_V = read(f, Float64)
                    read(f, Float64)#redundant sum_sq_B

                    N_corr = read(f, Int64)
                    corrs = Float64[]

                    for i = 1:N_corr
                        push!(corrs, read(f, Float64))
                    end

                    n_samples = read(f, Int64)

                    for i = 1:N_corr
                        corrs[i] /= (n_samples - i + 1)
                    end
                    avg_sq_grad_V = sum_sq_grad_V / n_samples

                    if proposal == "EM" # (STANDARD MALA)
                        Dhat = 1 - (β^2 * N * dt) * sum(corrs) / d
                    elseif rule == "METROPOLIS" #(METROPOLIS HMC)
                        Dhat = 1 - (β^2 * N * dt) * (sum(corrs)-0.5*first(corrs)) / d
                    else # (BARKER HMC)
                        Dhat = 1 - (β^2 * N * dt) * (sum(corrs)-first(corrs)) / (2d)
                    end
                    push!(est_Ds, Dhat)
                    push!(est_grad_V,first(corrs))
                    dt_range = dt * (0:(N_corr-1))
                   # plot(dt_range, corrs, label="")
                    #savefig("plots/autocorrelations/$(rule)_$(proposal)_$(lg_dt)_$(samp_count).png")

                end
                close(f)
                #println("GK estimates: ", est_Ds)
                push!(D_gk, last(est_Ds))
                push!(grad_V,last(est_grad_V))
            end

            if isfile("sd_increments_$(rule)_$(proposal)_$(lg_dt)")
                sd_file = true
                push!(dts_einstein, 10^lg_dt)
                dt = 10^lg_dt
                N_corr = ceil(Int64, T_corr / dt)
                T_inc = 10 * dt * N_corr
                f = open("sd_increments_$(rule)_$(proposal)_$(lg_dt)", "r")
                read(f,Float64) # discard first sample, systematically large
                incs = Float64[]

                while (!eof(f))
                    push!(incs, read(f, Float64))
                end

                close(f)
                #println("SD increments: ", incs)
                Dhat = sum(incs) / (length(incs) * 2 * d * N * T_inc)
                (rule == "BARKER") && (Dhat*=2) #adjust for rejection rate
                push!(D_einstein, Dhat)

            end

 
        end
    scatter!(normal_plot,dts_gk, D_gk, markershape=:xcross, label="GK $(rule) $(proposal)")
    scatter!(normal_plot,dts_einstein, D_einstein, markershape=:xcross, label="Einstein $(rule) $(proposal)")
    scatter!(eq_plot,dts_gk,grad_V,markershape=:xcross,label="$(rule) $(proposal)")
    #(any_file) && savefig("plots/diffusions_$(rule)_$(proposal).pdf")

    scatter!(log_plot,dts_gk, D_gk, markershape=:xcross, label="GK $(rule) $(proposal)", xlabel="Δt", ylabel="estimated diffusion")
    scatter!(log_plot,dts_einstein, D_einstein, markershape=:xcross, label="Einstein $(rule) $(proposal)")
    #(any_file) && savefig("plots/diffusions_$(rule)_$(proposal)_log.pdf")

    end

end

savefig(normal_plot,"plots/diffusions.pdf")
savefig(log_plot,"plots/diffusions_log.pdf")
savefig(eq_plot,"plots/eq_plot.pdf")