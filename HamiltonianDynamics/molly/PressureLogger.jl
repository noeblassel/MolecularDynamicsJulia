using Molly

struct PressureLogger{T}
    n_steps::Int
    pressures::Vector{T}
end

PressureLogger(T, n_steps::Integer)=PressureLogger(n_steps,T[])
PressureLogger(n_steps::Integer)=PressureLogger(Float32, n_steps)

function Base.show(io::IO, pl::PressureLogger)
    print(io, "PressureLogger{", eltype(eltype(pl.pressures)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.pressures), " pressures recorded")
end

function Molly.log_property!(logger::PressureLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        
    end
end