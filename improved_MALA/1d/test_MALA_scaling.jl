using Random, LinearAlgebra, Plots

include("integrators.jl")

Z_periodic=1.266066


potential(q::Float64) = sin(2π * q)
grad_potential(q::Float64) = 2π * cos(2π * q)
bc(q::Float64) = mod1(q,1.0)

"""potential(q::Float64)=0.5q^2
grad_potential(q::Float64)=q
bc(q::Float64)=q"""

M = 1000
q = rand(M)
qlims = (0.0,1.0)
N_bins = 200
q_range = range(qlims..., N_bins)
dq = q_range.step.hi

dt_ref = 0.01
N_max_iterations = 100_000

log_dts = range(-3, -1, 10)
dts = 10 .^ log_dts
#dts=[0.01]
Rs = zero(dts)
t_fin = N_max_iterations * minimum(dts)
q_hists = [zeros(N_bins) for dt in dts]
N_eq_steps = round(Int64, t_fin / dt_ref)
sim_eq = MALA(dt=dt_ref, T=1.0)
simulate!(q, sim_eq, N_eq_steps, potential, grad_potential, bc)#equilibriate

N_samples = 10
for i = 1:N_samples
    println("iteration $i/$N_samples")
    for (i, dt) in enumerate(dts)
        #N_steps = round(Int64, t_fin / dt)
        sim = MALA(dt=dt, T=1.0)
        simulate!(q, sim, N_max_iterations, potential, grad_potential, bc; record_hist=true, q_hist=q_hists[i], qbounds=qlims)
       # Rs[i] += 1 - sim.n_accepted / sim.n_total
    end
end
println(Rs / N_samples)
println(dts)

for (i,dt) in enumerate(dts)
    plot!(q_range,q_hists[i]/(dq*sum(q_hists[i])),linestyle=:dot,label="Δt=$(round(dt,digits=3))")
end
plot!(q->exp(-potential(q))/Z_periodic,qlims...,label="")
savefig("hists_mala_per.pdf")