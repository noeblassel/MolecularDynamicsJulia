using Molly


barker(p) = p / (1 + p)
sq_norm(v::Vector{SVector{3,Float64}}) = dot(v, v)

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
    β = ustrip(inv(T))
    MALA(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function Molly.simulate!(sys::System{D}, sim::MALA, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    grad_V = -forces(sys, neighbors; parallel=parallel)#compute the gradient, assumes all atoms have mass 1
    V = potential_energy(sys, neighbors)
    σ = sqrt(2 * sim.dt)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)
        dW = SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
        ## propose candidate position according to Euler-Maruyama scheme ##
        coords_tilde = sys.coords - sim.β * sim.dt * grad_V + σ * dW
        coords_tilde = wrap_coords_vec.(coords_tilde, (sys.box_size,))#enforce boundary conditions
        ## temporarily swap candidate position in the system's field to compute candidate potential and gradient ##
        sys.coords, coords_tilde = coords_tilde, sys.coords
        neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
        V_tilde = potential_energy(sys, neighbors)#candidate potential
        grad_V_tilde = -forces(sys, neighbors; parallel=parallel)#candidate gradient

        ## compute minus log of metropolis ratio ##
        α = sim.β * (V_tilde - V) # potential part
        +sq_norm(dW - sim.β * σ * (grad_V + grad_V_tilde) / 2) / 2 #reverse transition term
        -sq_norm(dW) / 2 # forward transition term

        ## accept / reject step ##
        U = rand(sim.rng)
        accepted = (sim.is_metropolis) ? (log(U) < -α) : (U < barker(exp(-α))) #Metropolis-Hastings  / Barker rule
        sim.n_total += 1

        if accepted
            sim.n_accepted += 1
            grad_V, grad_V_tilde = grad_V_tilde, grad_V #reuse candidate gradient as new gradient
            V = V_tilde #reuse candidate potential as new potential
        else
            sys.coords ,coords_tilde=coords_tilde,sys.coords #swap back original coordinates
        end
        neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    end

end

mutable struct MALA_HMC
    dt::Real
    β::Real

    is_metropolis::Bool

    rseed::UInt32
    rng::AbstractRNG

    n_accepted::Int64
    n_total::Int64
end

function MALA_HMC(; dt, T, is_metropolis=true, rseed=UInt32(round(time())), rng=MersenneTwister(rseed))
    β = ustrip(inv(T))
    MALA_HMC(dt, β, is_metropolis, rseed, rng, 0, 0)
end

function Molly.simulate!(sys::System{D}, sim::MALA_HMC, n_steps::Integer; parallel::Bool=true) where {D}
    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)

    V = potential_energy(sys, neighbors)

    h = sqrt(2 * sim.β * sim.dt) #effective timestep of the symplectic scheme
    σ = inv(sqrt(sim.β)) #noise scale

    sys.velocities=σ*SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D)))) #sample random momentum
    coords_tmp=zero(sys.coords)

    for i = 1:n_steps
        run_loggers!(sys, neighbors, i)
        

        # keep track of original position and Hamiltonian
        coords_tmp .= sys.coords 
        H=V+Molly.kinetic_energy_noconvert(sys)
        
        # position Verlet evolution
        sys.coords += h * sys.velocities / 2 #A
        sys.velocities += h * accelerations(sys,neighbors) #B
        sys.coords += h * sys.velocities / 2 #A

        sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))


        #compute new potential
        V_tilde=potential_energy(sys,neighbors)
        α=sim.β*(V_tilde+Molly.kinetic_energy_noconvert(sys)-H)

        ## accept / reject step ##
        U = rand(sim.rng)
        
        accepted = (sim.is_metropolis) ? (log(U) < -α) : (exp_minus_α = exp(-α);(U < (exp_minus_α / (1 + exp_minus_α))))
        sim.n_total += 1

        if accepted
            sim.n_accepted += 1
            neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel) #recompute neighbors
            V = V_tilde #reuse candidate potential as new potential
        else
            sys.coords .= coords_tmp
        end

        #sample new momentum for next iteration
        sys.velocities=σ*SVector{D}.(eachrow(randn(sim.rng, Float64, (length(sys), D))))
    end

end