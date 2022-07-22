using LinearAlgebra, StaticArrays
begin

    include("ToySystem.jl")
    include("ReflectingBoundaryCondition.jl")
    include("BAOAB.jl")
    include("GenParRep.jl")

    #domain, boundary condition and states definitions
    include("domain.jl")

    output_file="/libre/blasseln/MolecularDynamicsJulia/semaine_roscoff/transition_samples_gpr.out"
    gr_history_output_file="/libre/blasseln/MolecularDynamicsJulia/semaine_roscoff/gr_histories_gpr.out"

    dt=0.001

    #define methods for GPR algorithm
    function O(sys)
        i=get_state(sys)
        x,y=divrem(i,W)
        C=SVector(x*(L+2R),y*(L+2R))
        d=norm(C-sys.q)
        if d<=R
            return d
        else
            return 0.0
        end
    end

    function branch_replica(sys,rep_sys=nothing)
        if rep_sys == nothing
        return ToySystem(sys.q,sys.p,sys.V,sys.∇V,sys.boundary_condition!,[O])
        else #copy averages
            return ToySystem(sys.q,sys.p,sys.last_q,sys.V,sys.∇V,sys.boundary_condition!,[O],rep_sys.O_sums,rep_sys.sq_O_sums,rep_sys.clock)
        end
    end

    spawn_replica(sys)=ToySystem(sys.q,sys.p,sys.V,sys.∇V,sys.boundary_condition!)

    function get_gr_obs(sys)
        return sys.O_sums,sys.sq_O_sums
    end

    get_clock(sys)=sys.clock

    reset_clock!(sys)=(sys.clock=0)

    function output_transition(state,next_state,n_steps,gr_history=Vector{Vector{Float64}}[])
        t=n_steps*dt
        f=open(output_file,"a")
        println(f,"$state $next_state $t")
        close(f)
        println("transition from $state to $next_state")
        if length(gr_history)>0
            f=open(gr_history_output_file,"a")
            println(f,join(map(first,gr_history)," "))
            println(f,"*")
            close(f)
        end
    end

    algo=GenParRepAlgorithm(64,10,100,1,0.01,spawn_replica,branch_replica,get_gr_obs,get_state,get_clock,reset_clock!,output_transition)
    #
    p0=SVector(0.0,0.0)
    q0=rand(centers)

    p0=SVector(0.0,0.0)
    q0=SVector(0.0,0.0)

    obsA(sys)=norm(sys.q-C_A)
    obsB(sys)=norm(sys.q-C_B)

    V(q)=0.0
    ∇V(q)=[0.0,0.0]

    sys=ToySystem(q0,p0,V,∇V,bc!)
    sim=BAOABIntegrator(dt,1.0,4.0)

    f=open(output_file,"w")
    println(f,"state next_state transition_time")
    close(f)

    f=open(gr_history_output_file,"w")
    println("*")
    close(f)

    sample_transitions!(sys,algo,sim,100000)
end

