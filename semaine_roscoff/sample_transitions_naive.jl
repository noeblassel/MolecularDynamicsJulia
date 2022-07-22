begin
    using LinearAlgebra, StaticArrays

    include("ToySystem.jl")
    include("ReflectingBoundaryCondition.jl")
    include("BAOAB.jl")
    #domain, boundary condition and states definitions
    include("domain.jl")

    p0=SVector(0.0,0.0)
    q0=rand(centers)

    p0=SVector(0.0,0.0)
    q0=SVector(0.0,0.0)

    obsA(sys)=norm(sys.q-C_A)
    obsB(sys)=norm(sys.q-C_B)

    V(q)=0.0
    ∇V(q)=[0.0,0.0]

    sys=ToySystem(q0,p0,V,∇V,bc!,[])
    sim=BAOABIntegrator(0.005,1.0,0.7)

    check_every=100
    last_transition_time=0.0
    output_file="/libre/blasseln/MolecularDynamicsJulia/semaine_roscoff/transition_samples_naive.out"
    f=open(output_file,"w")
    println(f,"state next_state transition_time")
    close(f)

    n_transitions=0
    while n_transitions<100000
        last_state=get_state(sys)
        simulate!(sys,sim,check_every)
        state=get_state(sys)
        if state!=last_state
            println("transition from state $last_state to state $state")
            g=open(output_file,"a")
            t=sim.dt*sys.clock
            println(g,"$last_state $state $t")
            close(g)
            sys.clock=0
            global n_transitions+=1
        end
    end
end

