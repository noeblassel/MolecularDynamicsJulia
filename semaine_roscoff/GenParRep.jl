using Base.Threads, Statistics

struct GenParRepAlgorithm{SR,BR,GR,CS,GT,RC,OT}
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
    get_clock::GT #method to get the simulation time of a given system
    reset_clock!::RC #method to reset the clock of a system
    output_transition::OT #method to output a transition event, of the form (old_state, next_state, transition_time)

end

function sample_transitions!(master_system,algorithm::GenParRepAlgorithm,simulator,num_transitions::Int)
    gr_checks_per_state_checks= div(algorithm.n_steps_check,algorithm.n_steps_gr)
    n_transitions = 0
    current_master_state = algorithm.get_state(master_system)

    while n_transitions<num_transitions
        algorithm.reset_clock!(master_system)
        slave_systems=[algorithm.branch_replica(master_system) for i=1:algorithm.n_replicas] #spawn replicas for Flemming-Viot equilibriation

        gr_history=Vector{Float64}[]
        transitioned = false
        equilibriated = false

        while !transitioned
            simulate!(master_system,simulator,algorithm.n_steps_check)
            master_state = algorithm.get_state(master_system)

            if master_state != current_master_state #trigger transition output and break
                #println("Master particle has exited")
                n_transitions+=1
                algorithm.output_transition(master_state,current_master_state,algorithm.get_clock(master_system))
                transitioned=true
                current_master_state=master_state
            elseif !equilibriated #Flemming-Viot

                for substep=1:gr_checks_per_state_checks
                    @threads for replica ∈ slave_systems #Multithread Flemming-Viot, later should be upgraded to a more sophisticated architecture
                        simulate!(replica,simulator,algorithm.n_steps_gr)
                    end
                end

                dead_replicas_ix = Int[]

                Os=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]
                O2s=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]

                for (i,replica) ∈ enumerate(slave_systems)
                    replica_state=algorithm.get_state(replica)
                    
                    if replica_state != master_state #kill replica
                        #println("Replica $i is killed after $n_clock steps")
                        push!(dead_replicas_ix,i)
                    end
                end
                
                if length(dead_replicas_ix)>0 #branch if dead replicas
                    alive_replicas_ix=setdiff(1:algorithm.n_replicas,dead_replicas_ix)
                    for i ∈ dead_replicas_ix
                        new_ix=rand(alive_replicas_ix)
                        slave_systems[i]=algorithm.branch_replica(slave_systems[new_ix],slave_systems[i])
                    end
                end

                Os=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]
                sq_Os=[zeros(algorithm.n_gr_observables) for i=1:algorithm.n_replicas]

                for (i,replica) ∈ enumerate(slave_systems)
                    O,sq_O=algorithm.get_gr_obs(replica)
                    n_clock=algorithm.get_clock(replica)
                    Os[i] = O/n_clock
                    sq_Os[i] = sq_O/n_clock
                end
                #println("Obar: $(first(Obar)), O2bar: $(first(O2bar))")

                Obar = mean(Os)
                num=mean(sq_Os[i]-2Os[i] .* Obar + Obar .^ 2 for i=1:algorithm.n_replicas)
                denom=mean(sq_Os[i] - Os[i] .^ 2 for i=1:algorithm.n_replicas)
                GR = num ./ denom
                push!(gr_history,GR)
                #println(last(gr_history))
                #println("num: $num, denom: $denom, GR: $GR")

                if all( gr < 1+algorithm.gr_tol for gr ∈ GR)
                    equilibriated = true
                end

            else #Parallel replicas
                slave_systems=[algorithm.spawn_replica(replica) for replica ∈ slave_systems]
                @threads for i=1:algorithm.n_replicas
                    replica=slave_systems[i]
                    simulate!(replica,simulator,algorithm.n_steps_check)
                    replica_state=algorithm.get_state(replica)
                    if replica_state != current_master_state
                        clock_replica= algorithm.get_clock(replica)
                        n_transition = (algorithm.n_replicas-1)*clock_replica + i
                        algorithm.output_transition(current_master_state,replica_state,n_transition,gr_history)
                        #println(gr_history)
                        master_system=replica
                        transitioned=true
                        n_transitions+=1
                        current_master_state=replica_state
                    end
                end
            end
        end
    end
end
