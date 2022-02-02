using Molly

struct PressureLogger{T}
    n_steps::Int
    pressures::Vector{Vector{T}}
end

function PressureLogger(T, n_steps::Integer; dims::Integer=3)
    return PressureLogger(n_steps,
                            Array{SArray{Tuple{dims}, T, 1, dims}, 1}[])
end

function PressureLogger(n_steps::Integer; dims::Integer=3)
    @assert dims==3 "pressure logger only available for 3D systems"
    return PressureLogger(typeof(one(Float32)u"Pa"), n_steps; dims=dims)
end

function Base.show(io::IO, pl::PressureLogger)
    print(io, "PressureLogger{", eltype(eltype(pl.pressures)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.pressures), " pressures recorded")
end

function log_property!(logger::PressureLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        
    end
end