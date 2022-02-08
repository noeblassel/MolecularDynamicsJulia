struct HamiltonianLogger{T}
    n_steps::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer)=HamiltonianLogger(n_steps,T[])
HamiltonianLogger(n_steps::Integer)=HamiltonianLogger(Float32, n_steps)

function Base.show(io::IO, hl::HamiltonianLogger)
    print(io, "HamiltonianLogger{", eltype(eltype(pl.energies)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.energies), " energies recorded")
end

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies,Molly.kinetic_energy_noconvert(s)+Molly.potential_energy(s))
    end
end