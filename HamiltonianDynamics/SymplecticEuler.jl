
###from Molly source, Molly.jl/src/simulators.jl
function remove_molar(x)
    fx = first(x)
    if dimension(fx) == u"𝐋 * 𝐍^-1 * 𝐓^-2"
        T = typeof(ustrip(fx))
        return x / T(Unitful.Na)
    else
        return x
    end
end
###


struct SymplecticEulerA{T,C}
    dt::T
    coupling::C
end

SymplecticEulerA(; dt, coupling=NoCoupling()) = SymplecticEulerA(dt, coupling)

struct SymplecticEulerB{T,C}
    dt::T
    coupling::C
end

SymplecticEulerA(; dt, coupling=NoCoupling()) = SymplecticEulerA(dt, coupling)


function simulate!(sys::System{D, true},
    sim::SymplecticEulerA,
    n_steps::Integer;
    parallel::Bool=true) where {D, S}

@assert sys.coupling <:Molly.NoCoupling "SymplecticEulerA does not currently support coupling"

if any(inter -> !inter.nl_only, values(sys.general_inters))
neighbors_all = all_neighbors(length(sys))
else
neighbors_all = nothing
end
neighbors = find_neighbors(sys, sys.neighbor_finder)

for step_n in 1:n_steps
run_loggers!(sys, neighbors, step_n)

sys.coords += sys.velocities .* sim.dt
sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))
sys.velocities += remove_molar.(accelerations(sys, sys.coords, sys.atoms, neighbors, neighbors_all)) .* sim.dt

if step_n != n_steps
neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n)
end
end
return sys
end


function simulate!(sys::System{D, true},
    sim::SymplecticEulerB,
    n_steps::Integer;
    parallel::Bool=true) where {D, S}

@assert sys.coupling <:Molly.NoCoupling "SymplecticEulerB does not currently support coupling"

if any(inter -> !inter.nl_only, values(sys.general_inters))
neighbors_all = all_neighbors(length(sys))
else
neighbors_all = nothing
end
neighbors = find_neighbors(sys, sys.neighbor_finder)

for step_n in 1:n_steps
run_loggers!(sys, neighbors, step_n)

sys.velocities += remove_molar.(accelerations(sys,sys.coords,sys.atoms,neighbors,neighbors_all))
sys.coords+=sys.velocities .* sim.dt
sys.coords = wrap_coords_vec.(sys.coords, (sys.box_size,))

if step_n != n_steps
neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n)
end
end
return sys
end