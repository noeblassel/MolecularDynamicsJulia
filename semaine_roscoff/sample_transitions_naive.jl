begin
    using LinearAlgebra, StaticArrays

    include("ToySystem.jl")
    include("ReflectingBoundaryCondition.jl")
    include("BAOAB.jl")
    #domain, boundary condition and states definitions
    include("domain.jl")

    p0=SVector(0.0,0.0)
    q0=rand(centers)

    
    obsA(sys)=norm(sys.q-C_A)
    obsB(sys)=norm(sys.q-C_B)
    
    V(q)=0.0
    ∇V(q)=[0.0,0.0]
    
    sys=ToySystem(q0,p0,V,∇V,bc!,[])
    current_state=get_state(sys)
    sim=BAOABIntegrator(0.001,1.0,4.0)

    check_every=600
    last_transition_time=0.0
    output_file="/libre/blasseln/MolecularDynamicsJulia/semaine_roscoff/transition_samples_naive.out"
    # output_file="/home/noeblassel/Documents/stage_CERMICS_2022/semaine_roscoff/transition_samples_naive.out"
    f=open(output_file,"w")
    println(f,"state next_state transition_time")
    close(f)

    n_transitions=0
    clock=0
    while n_transitions<100000
        simulate!(sys,sim,check_every)
        clock+=check_every
        state=get_state(sys,current_state)
        if state!=current_state
            println("transition from state $current_state to state $state")
            g=open(output_file,"a")
            t=sim.dt*clock
            println(g,"$current_state $state $t")
            close(g)
            global clock=0
            global n_transitions+=1
            global current_state=state
        end
    end
end

