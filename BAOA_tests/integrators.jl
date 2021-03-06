abstract type BoundaryCondition end

struct PeriodicBoundaryCondition <: BoundaryCondition
    L::Float64
end

struct InfiniteBox <: BoundaryCondition
    L::Float64
end

InfiniteBox() = InfiniteBox(Inf)

struct LangevinSplitting
    dt::Float64
    γ::Float64
    β::Float64

    rseed::UInt32
    rng::AbstractRNG

    splitting::AbstractString
    bc::BoundaryCondition
end

function LangevinSplitting(; dt, γ, T, splitting, rseed=UInt32(round(time())), rng=MersenneTwister(rseed), bc=InfiniteBox())
    β = inv(T)
    @assert (all(x ∈ "ABO" for x ∈ splitting) && all(x ∈ splitting for x ∈ "ABO")) "Invalid splitting descriptor: use only and all letters A, B and O."
    LangevinSplitting(dt, γ, β, rseed, rng, splitting, bc)
end


function simulate!(p_vec::Vector{Float64}, q_vec::Vector{Float64}, ∇V::Function, ∇K::Function, hist::Array{Int64,2}, qlims::Tuple{Float64,Float64}, plims::Tuple{Float64,Float64}, sim::LangevinSplitting, n_steps::Integer)

    α_eff = exp(-sim.γ * sim.dt / count('O', sim.splitting))
    σ_eff = sqrt((1 - α_eff^2) / sim.β)
    grad_vec = ∇V.(q_vec, (sim.bc.L,))
    effective_dts = [sim.dt / count(c, sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    steps = []
    arguments = []

    for (j, op) in enumerate(sim.splitting)
        if op == 'A'
            push!(steps, A_step!)
            push!(arguments, (q_vec, p_vec, effective_dts[j], ∇K, sim.bc))
        elseif op == 'B'
            push!(steps, B_step!)
            push!(arguments, (q_vec, p_vec, effective_dts[j], grad_vec, ∇V, force_computation_steps[j], sim.bc))
        elseif op == 'O'
            push!(steps, O_step!)
            push!(arguments, (p_vec, α_eff, σ_eff, sim.rng))
        end
    end

    step_arg_pairs = zip(steps, arguments)

    for step_n = 1:n_steps
        update_hist2D!.((hist,), q_vec, p_vec, (qlims,), (plims,))
        for (step!, args) = step_arg_pairs
            step!(args...)
        end

        (step_n % 100000 == 0) && (println(step_n, "/", n_steps, " steps done."); flush(stdout))
    end
end

function simulate2D!(p_vec::Vector{Float64}, q_vec::Vector{Float64}, force::Function, hist::Array{Int64,2}, M::Vector{Float64}, plims::Tuple{Float64,Float64}, sim::LangevinSplitting, n_steps::Integer)
    l = Int64(size(q_vec, 1) // size(M, 1))
    M_full = repeat(M, l)
    M_inv = inv.(repeat(M, l))
    # println(M_full)
    α_eff = zero(q_vec)
    σ_eff = zero(q_vec)
    @. α_eff = exp(-sim.γ * sim.dt * M_inv / count('O', sim.splitting))
    @. σ_eff = sqrt((1 - α_eff^2) * M_full / sim.β)
    # println(α_eff)
    # println(σ_eff)
    force_vec = reduce(vcat, force.(q_vec[1:2:end], q_vec[2:2:end], (sim.bc.L,)))
    # println(force_vec)
    effective_dts = [sim.dt / count(c, sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    steps = []
    arguments = []

    for (j, op) in enumerate(sim.splitting)
        if op == 'A'
            push!(steps, A_step2D!)
            push!(arguments, (q_vec, p_vec, M_inv, effective_dts[j], sim.bc))
        elseif op == 'B'
            push!(steps, B_step2D!)
            push!(arguments, (q_vec, p_vec, effective_dts[j], force_vec, force, force_computation_steps[j], sim.bc))
        elseif op == 'O'
            push!(steps, O_step_vec!)
            push!(arguments, (p_vec, α_eff, σ_eff, sim.rng))
        end
    end

    step_arg_pairs = zip(steps, arguments)
    #println("---------------")
    for step_n = 1:n_steps
        update_hist2D!.((hist,), p_vec[1:2:end], p_vec[2:2:end], (plims,), (plims,))
        for (step!, args) = step_arg_pairs
            step!(args...)
        end
        # println("p_vec:")
        # println(p_vec)
        # println("q_vec:")
        # println(q_vec,"\n-----------")

        (step_n % 100000 == 0) && (println(step_n, "/", n_steps, " steps done."); flush(stdout))
    end
end

function O_step!(p_vec::Vector{Float64}, α_eff::Float64, σ_eff::Float64, rng::AbstractRNG)
    p_vec .= α_eff * p_vec + σ_eff * randn(rng, Float64, size(p_vec))
end

function O_step_vec!(p_vec::Vector{Float64}, α_eff::Vector{Float64}, σ_eff::Vector{Float64}, rng::AbstractRNG)
    p_vec .= α_eff .* p_vec + σ_eff .* randn(rng, Float64, size(p_vec))
end

function A_step!(q_vec::Vector{Float64}, p_vec::Vector{Float64}, dt_eff::Float64,∇K::Function, bc::BoundaryCondition=InfiniteBox())
    q_vec .+= ∇K.(p_vec) * dt_eff
    (isa(bc, PeriodicBoundaryCondition)) && (q_vec .= mod1.(q_vec, (bc.L,)))
end

function A_step2D!(q_vec::Vector{Float64}, p_vec::Vector{Float64}, M_inv::Vector{Float64}, dt_eff::Float64, bc::BoundaryCondition=InfiniteBox())
    q_vec .+= (M_inv .* p_vec) * dt_eff
    (isa(bc, PeriodicBoundaryCondition)) && (q_vec .= mod1.(q_vec, (bc.L,)))
end

function B_step!(q_vec::Vector{Float64}, p_vec::Vector{Float64}, dt_eff::Float64, grad_vec::Vector{Float64}, ∇V::Function, compute_forces::Bool, bc::BoundaryCondition=InfiniteBox())
    compute_forces && (grad_vec .= ∇V.(q_vec, (bc.L,)))
    p_vec .-= dt_eff * grad_vec
end

function B_step2D!(q_vec::Vector{Float64}, p_vec::Vector{Float64}, dt_eff::Float64, force_vec::Vector{Float64}, force_func::Function, compute_forces::Bool, bc::BoundaryCondition=InfiniteBox())
    compute_forces && (force_vec .= reduce(vcat, force.(q_vec[1:2:end], q_vec[2:2:end], (sim.bc.L,))))
    p_vec .+= dt_eff * force_vec
end