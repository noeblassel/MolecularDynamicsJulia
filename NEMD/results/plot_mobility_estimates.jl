using Plots,LinearAlgebra
#TODO add error bars using block averaging

path_orig="/libre/blasseln/MolecularDynamicsJulia/NEMD/*.out"
path_end="."

run(`scp clustern15 $path_orig $path_end`)
run(`scp clustern16 $path_orig $path_end`)

ηs=0.1:0.1:1.0

methods=["SINGLE","COLOR"]

for m in methods
    println(m)
       Rs=[]
       for η in ηs
        println("\t",η)
        f=open("mobility_estimates$(η)_$(m).out","r")
        n_samps=0
        sum_R=0.0
        while !eof(f)
          sum_R+=read(f,Float64)
          n_samps+=1
        end
        close(f)
        push!(Rs,sum_R/n_samps)
       end
  (m=="COLOR") && (Rs*=sqrt(1000))
  a=inv(dot(ηs,ηs))*dot(ηs,Rs)#least squares slope fit
  println(a)
  scatter(ηs,Rs,markershape=:xcross,label=m,color=:blue,legend=:topleft)
  plot!(x->a*x,0,last(ηs),linestyle=:dot,color=:red,label="fit")
  savefig("$(m).pdf")
end
