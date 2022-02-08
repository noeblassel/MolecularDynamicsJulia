struct KineticEnergyLoggerNoDims{T}
    n_steps::Int
    energies::Vector{T}
end

KineticEnergyLoggerNoDims(T, n_steps::Integer)=KineticEnergyLoggerNoDims(n_steps,T[])
KineticEnergyLoggerNoDims(n_steps::Integer)=KineticEnergyLoggerNoDims(Float32, n_steps)

function Base.show(io::IO, pl::KineticEnergyLoggerNoDims)
    print(io, "KineticEnergyLoggerNoDims{", eltype(eltype(pl.energies)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.energies), " energies recorded")
end

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies,Molly.kinetic_energy_noconvert(s))
    end
end