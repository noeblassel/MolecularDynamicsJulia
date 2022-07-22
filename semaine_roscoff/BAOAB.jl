struct BAOABIntegrator
    dt::Float64
    γ::Float64
    β::Float64
end

function simulate!(sys::ToySystem{D},sim::BAOABIntegrator,n_steps::Int) where {D}

    grad_V=sys.∇V(sys.q)

    α=exp(-sim.γ * sim.dt)
    σ=sqrt((1-α^2)/sim.β)

    for i=1:n_steps
        sys.last_q=sys.q
        
        if length(sys.observables) != 0
            observable=Float64[obs(sys) for obs in sys.observables]

            sys.O_sums += observable
            sys.sq_O_sums += observable .^ 2
        end

        sys.p -= sim.dt*sys.∇V(sys.q)/2
        sys.q += sim.dt*sys.p/2

        sys.boundary_condition!(sys)

        grad_V=sys.∇V(sys.q)
        sys.p = α*sys.p + σ*randn(D)
        sys.last_q=sys.q
        sys.q+= sim.dt*sys.p/2
        sys.boundary_condition!(sys)
        sys.p -= sim.dt*sys.∇V(sys.q)/2
    end
end