barker(p::Float64) = p / (1 + p)

mutable struct MALA
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
end

function MALA(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = inv(T)
    MALA(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function simulate!(q::Vector{Float64}, sim::MALA, n_steps::Integer, potential::F_V, grad_potential::F_DV, bc::F_BC; record_hist::Bool=false, q_hist::Vector{Float64}=Float64[], qbounds::Tuple{Float64,Float64}=(0.0, 0.0)) where {F_V,F_DV,F_BC}
    M = length(q) #number of parallel replicas
    grad_V = grad_potential.(q)
    V = potential.(q)
    σ = sqrt(2 * sim.dt)
    for i = 1:n_steps
        record_hist && update_hist!.((q_hist,),(qbounds,), q)
        G = randn(sim.rng, Float64, M)
        qtilde = q - sim.β * sim.dt * grad_V + σ * G
        displacement = qtilde - q
        qtilde = bc.(qtilde)
        Vtilde = potential.(qtilde)
        grad_Vtilde = grad_potential.(qtilde)
        alpha = sim.β * (Vtilde - V) + ((sim.β * sim.dt * grad_Vtilde - displacement) .^ 2) / (4 * sim.dt) - ((displacement + sim.β * sim.dt * grad_V) .^ 2) / (4 * sim.dt)
        U = rand(sim.rng, Float64, M)
        accepted = log.(U) .< -alpha

        q[accepted] .= qtilde[accepted]
        V[accepted] .= Vtilde[accepted]
        grad_V[accepted] .= grad_Vtilde[accepted]
        sim.n_accepted += sum(accepted)
        sim.n_total+=M
    end
end


mutable struct EM
    dt::Real
    β::Real

    rseed::UInt32
    rng::AbstractRNG
end

function EM(; dt, T, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = inv(T)
    EM(dt, β, rseed, rng)
end

function simulate!(q::Vector{Float64}, sim::EM, n_steps::Integer, potential::F_V, grad_potential::F_DV, bc::F_BC; record_hist::Bool=false, q_hist::Vector{Float64}=Float64[], qbounds::Tuple{Float64,Float64}=(0.0, 0.0)) where {F_V,F_DV,F_BC}
    M = length(q) #number of parallel replicas
    grad_V = grad_potential.(q)
    V = potential.(q)
    σ = sqrt(2 * sim.dt)
    #println("-----------------------------")
    for i = 1:n_steps
        record_hist && update_hist!.((q_hist,),(qbounds,), q)
        G = randn(sim.rng, Float64, M)
        q = q - sim.β * sim.dt * grad_V + σ * G

    end
end

function update_hist!(qhist::Vector{Float64}, qbounds::Tuple{Float64,Float64}, q::Float64)
    ix=1+floor(Int64,length(qhist)*(q-first(qbounds))/(last(qbounds)-first(qbounds)))
    (0<ix) && (length(qhist)>=ix) && (qhist[ix]+=1)
end