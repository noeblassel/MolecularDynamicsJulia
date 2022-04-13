#!/bin/env julia


using Plots

cmd="./fetch_mobility_estimates.sh"
arg=lowercase(ARGS[1])

n_regr_points=3

function read_last(file)#stolen from stack overflow
    open(file) do io
      seekend(io)
      seek(io, position(io) - 2)
      while Char(peek(io)) != '\n'
        seek(io, position(io) - 1)
      end
      read(io, Char)
      return read(io, String)
    end
  end

run(`$cmd $arg`)

data_folder="./mobility_computations_$(arg)/"
data_files=readdir(data_folder)

ηs=Float64[]
Rs=Float64[]

for file in data_files
    m=match(r"mobility_estimates(.+)\.out",file)
    push!(ηs,parse(Float64,first(m.captures)))
    l=read_last(data_folder*file)
    (_,ρ)=split(l)
    push!(Rs,parse(Float64,ρ))
end

normalization= (arg=="color")  ? 1000.0 : 1.0
Rs/=normalization

scatter(ηs,Rs,markershape=:xcross,label="$(arg) drift", xlabel="η",ylabel="R(η)",legend=:top,color=:blue)

perm=sortperm(ηs)
ηs=ηs[perm]
Rs=Rs[perm]

#fit R=aη
a=sum(Rs[1:n_regr_points])*inv(sum(ηs[1:n_regr_points]))
f(x)=a*x

#mobility estimation
D= (arg == "color") ? (999/1000)*(a + 1/999) : a

println("Estimated mobility : $(D)")
plot!(f,0,maximum(ηs),linestyle=:dot,color=:red,label="y=$(a)x")

savefig("plot_$(arg).pdf")