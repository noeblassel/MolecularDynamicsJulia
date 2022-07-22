using Base.Threads, Statistics

struct GenParRepAlgorithm{SR,BR,GR,CS,OT}
    n_replicas::Int

    n_steps_gr::Int #simulation iterations between checks of GR statistics convergence
    n_steps_check::Int #simulation iterations between checks of system states
    n_gr_observables::Int #number of observables for Gelman-Rubin
    gr_tol::Float64 #tolerance for assessing convergence to the QSD via Gelman-Rubin statistic

    #functions
    spawn_replica::SR #method to spawn a replica from a master system
    branch_replica::BR #method to branch a replica for Flemming-Viot
    get_gr_obs::GR #method to get running sum and running sum of squares for Gelman-Rubin observables
    get_state::CS #method to check the state of the system
    output_transition::OT #method to output a transition event, of the form (old_state, next_state, transition_time)

end

function sample_transitions!(master_system,algorithm::GenParRepAlgorithm,simulator,num_transitions::Int)
    gr_checks_per_state_checks= div(algorithm.n_steps_check,algorithm.n_steps_gr)
    n_transitions = 0
    current_master_state=algorithm.get_state(master_system)
    
    while n_transitions<num_transitions
        #println("$n_transitions transitions observed so far, master is in state $current_master_state")
        flemming_viot_clock = 0
        parallel_exit_clock = 0
        slave_systems=[algorithm.branch_replica(master_system) for i=1:algorithm.n_replicas] #spawn replicas for Flemming-Viot equilibriation

        transitioned = false
        
        ### equilibriation step
        equilibriated = false
        gr_history=Vector{Float64}[]
        
        #println("Equilibriating to QSD")

        while !equilibriated
            simulate!(master_system,simulator,algorithm.n_steps_check)
            master_state = algorithm.get_state(master_system,current_master_state)

            if master_state != current_master_state
                n_transitions +=1
                #println("Master particle exits from $current_master_state to $master_state")
                algorithm.output_transition(current_master_state,master_state,flemming_viot_clock+algorithm.n_steps_check)
                transitioned = true
                current_master_state=master_state
                break
            end

            dead_replicas_ix = Int[]

            for (i,replica) ∈ enumerate(slave_systems)
                replica_state=algorithm.get_state(replica,current_master_state)
                
                if replica_state != master_state #Detect replica exits
                    push!(dead_replicas_ix,i)
                end
            end

            if length(dead_replicas_ix)>0 #Branching step to kill exited replicas
                #println("$(length(dead_replicas_ix)) have exited")
                alive_replicas_ix=setdiff(1:algorithm.n_replicas,dead_replicas_ix)
                for i ∈ dead_replicas_ix
                    new_ix=rand(alive_replicas_ix) #pick random active particle
                    slave_systems[i]=algorithm.branch_replica(slave_systems[new_ix],slave_systems[i])#replace dead replica by an independent copy of the picked particle
                end
            end

            for substep=1:gr_checks_per_state_checks #Check convergence to QSD via Gelman-Rubin statistic
                @threads for i=1:algorithm.n_replicas #@threads eventually becomes a more sophisticated architecture
                    simulate!(slave_systems[i],simulator,algorithm.n_steps_gr)
                end
                flemming_viot_clock += algorithm.n_steps_gr

                Os=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]
                sq_Os=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]

                for (i,replica) ∈ enumerate(slave_systems)
                    O_i,sq_O_i=algorithm.get_gr_obs(replica)
                    Os[i] = O_i/flemming_viot_clock
                    sq_Os[i] = sq_O_i/flemming_viot_clock
                end

                Obar = mean(Os)
                num=mean(sq_Os[i]-2Os[i] .* Obar + Obar .^ 2 for i=1:algorithm.n_replicas)
                denom=mean(sq_Os[i] - Os[i] .^ 2 for i=1:algorithm.n_replicas)
                GR = (num ./ denom) .- 1.0
                push!(gr_history,GR)

                if all( gr < algorithm.gr_tol for gr ∈ GR)
                    equilibriated = true
                    break
                end

            end

        end
        #println("Equilibriated in $flemming_viot_clock steps, going to parallel step")
        ##println("gr history : ", map(first,gr_history))
        while !transitioned #Replicas have equilibriated
            slave_systems=[algorithm.spawn_replica(replica) for replica ∈ slave_systems]
            exited_ix=[Int[] for i=1:nthreads()]
            
            @threads for i=1:algorithm.n_replicas
                replica=slave_systems[i]
                simulate!(replica,simulator,algorithm.n_steps_check)
                replica_state=algorithm.get_state(replica,current_master_state)
                if replica_state != current_master_state
                    
                    push!(exited_ix[threadid()],i)
                    break
                end
            end

            parallel_exit_clock += algorithm.n_steps_check

            exited_ix=vcat(exited_ix...)

            if length(exited_ix) >0
                ix = rand(exited_ix)
                #println("Replica $ix has won the race in $parallel_exit_clock steps")
                n_exit_steps = flemming_viot_clock + (algorithm.n_replicas-1)*parallel_exit_clock + ix
                master_system=slave_systems[ix]
                next_state=get_state(master_system,current_master_state)
                algorithm.output_transition(current_master_state,next_state,n_exit_steps,gr_history)
                current_master_state=next_state
                transitioned=true
                n_transitions+=1
            end
        end #while !transitioned
    end #while n_transitions<num_transitions
    
end