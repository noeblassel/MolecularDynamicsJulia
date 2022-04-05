using LsqFit
using ForwardDiff

"Usage: julia fit_correlation.jl dt Y_file N_sinusoids"

dt=parse(Float64,ARGS[1])
Y=parse.(Float64,split(readline(ARGS[2]),","))
N_sinusoids=parse(Int64,ARGS[3])

X=(0:length(Y)-1)*dt
X=[X...]
model(x,p)=exp.(-p[1]*x).*sum(p[3i-1]*sin.(p[3i]*x .+ p[3i+1]) for i=1:N_sinusoids)
#parameters
P0=ones(3*N_sinusoids+1)

fit=curve_fit(model,X,Y,P0)

dump(fit)

f=open("fit_"*ARGS[2],"w")
println(f,model.(X,(coef(fit),)))
close(f)
