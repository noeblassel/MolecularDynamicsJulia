dt=1e-3
T=0.1
β=inv(T)

f=open("results/autocorrelation_history.out","r")
N=read(f,Int64)
sum_A=zeros(3N)

for i=1:3N
    sum_A[i]=read(f,Float64)
end

N=read(f,Int64)
sum_B=zeros(3N)

for i=1:3N
    sum_B[i]=read(f,Float64)
end

L=read(f,Int64)
C=zeros(L)

for i=1:L
    C[i]=read(f,Float64)
end

n_timesteps=read(f,Int64)
close(f)

for i=eachindex(C)
    C[i]/=(n_timesteps-i+1)
end
## C= [sum(V_iΔt. V_0)/N^2, i=1..L]
using Plots, LinearAlgebra

t_range=0:dt:(dt*(L-1))
plot(t_range,C,label="")
savefig("velocity_autocorrelation.pdf")

I=0
Int_curve=[]
for (c1,c2)=zip(C,C[2:end])
    global I+=(c1+c2)*dt/2
    push!(Int_curve,I)
end
println("Estimated transport coefficient :",β*I)

plot(t_range[1:2000],C[1:2000],label="")
savefig("short_time_autocorrelation.pdf")
plot(t_range[1:2000],Int_curve[1:2000])
savefig("int_curve.pdf")
