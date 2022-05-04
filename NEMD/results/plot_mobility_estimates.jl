using Plots,LinearAlgebra

ηs=[0.01,0.02,0.03,0.04,0.05]
methods=["SINGLE","COLOR"]
for m in methods
       Rs=[]
       for η in ηs
        f=open("mobility_estimates$(η)_$(m).out","r")
        read(f,Int64)
        push!(Rs,read(f,Float64))
        close(f)
       end

  a=inv(dot(ηs,ηs))*dot(ηs,Rs)#least squares slope fit
  println(a)
  scatter(ηs,Rs,markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(x->a*x,0,last(ηs),linestyle=:dot,color=:red,label="fit")
  savefig("$(m).pdf")
end
