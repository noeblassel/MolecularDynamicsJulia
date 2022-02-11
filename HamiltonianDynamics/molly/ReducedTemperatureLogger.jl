struct ReducedTemperatureLogger{T}
    n_steps::Int
    temperatures::Vector{T}
end

ReducedTemperatureLogger(T, n_steps::Integer)=ReducedTemperatureLogger(n_steps,T[])
ReducedTemperatureLogger(n_steps::Integer)=ReducedTemperatureLogger(Float64, n_steps)

function Base.show(io::IO, pl::ReducedTemperatureLogger)
    print(io, "ReducedTemperatureLogger{", eltype(eltype(pl.pressures)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.pressures), " temperatures recorded")
end

function Molly.log_property!(logger::ReducedTemperatureLogger, s::System, neighbors=nothing, step_n::Integer=0)
    N=length(s)
    if step_n % logger.n_steps == 0
        ke=Molly.kinetic_energy_noconvert(s)
        T=2*ke/(3*(N-1))
        push!(logger.temperatures,T)
    end
end