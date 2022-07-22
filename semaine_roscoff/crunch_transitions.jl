base_path="/home/noeblassel/Documents/stage_CERMICS_2022/semaine_roscoff/"
input_file="transition_samples_naive.out"

W=3
H=3
transitions=zeros(W*H,W*H)
transition_times=zeros(W*H,W*H)
exit_times=[Float64[] for i=1:W*H]
n_transitions=0

for l in readlines(joinpath(base_path,input_file))[2:end]
    state,next_state,tau=split(l," ")
    state=parse(Int,state)
    next_state=parse(Int,next_state)
    tau=parse(Float64,tau)

    global n_transitions+=1
    transitions[state+1,next_state+1]+=1.0
    transition_times[state+1,next_state+1]+=tau
    push!(exit_times[state+1],tau)
end

#normalize_quantities

n_transitions_from=Int[]

for i=1:W*H
    transitions_from_i=sum(transitions[i,:])
    push!(n_transitions_from,transitions_from_i)
    transitions[i,:] /= transitions_from_i
    transition_times[i,:] /= transitions_from_i
end

f=open(joinpath(base_path,"digested_"*input_file),"w")
println(f,"transition matrix:")
println(f,transitions)
println(f,"transition times matrix:")
println(f,transition_times)

for i=1:W*H
    println(f,"\tState: $(i-1)")
    println(f,"\tmean_transition_time:")
    println(f,"\t",sum(exit_times[i])/n_transitions_from[i])
    println(f,"\texit_times:")
    println(f,exit_times[i])
end

close(f)


