using Plots

dt=5e-3
N=1000
T=2.5
β=inv(T)
avg_ρ=0.0
t_between_samples=1000.0

for i=0:10
    r=parse.(Float64,readlines("msd_history$(i).out"))
    n_steps=length(r)
    t_fin=n_steps*t_between_samples
    plot!(t_between_samples:t_between_samples:t_fin,r,label="")
    global avg_ρ+=β*last(r)/(6*t_fin)
end

println("estimated ρ: $(avg_ρ/11)")

savefig("superpose_msd.pdf")