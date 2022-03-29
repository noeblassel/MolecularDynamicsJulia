###Total energy logger

struct HamiltonianLogger{T}
    n_steps::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer) = HamiltonianLogger(n_steps, T[])
HamiltonianLogger(n_steps::Integer) = HamiltonianLogger(Float32, n_steps)

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, Molly.kinetic_energy_noconvert(s) + Molly.potential_energy(s, neighbors))
    end
end

###Dimensionless kinetic energy logger--- Molly version (as of now) requires physical units

struct KineticEnergyLoggerNoDims{T}
    n_steps::Int
    energies::Vector{T}
end

KineticEnergyLoggerNoDims(T, n_steps::Integer) = KineticEnergyLoggerNoDims(n_steps, T[])
KineticEnergyLoggerNoDims(n_steps::Integer) = KineticEnergyLoggerNoDims(Float64, n_steps)

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, Molly.kinetic_energy_noconvert(s))
    end
end

###State logger (writes state of system to external file --- NOT GENERAL)

struct StateLogger
    n_steps::Int
    prefix::AbstractString
end

StateLogger(n_steps::Integer) = StateLogger(n_steps, "logfile")
StateLogger(n_steps::Integer, file_prefix::AbstractString) = StateLogger(n_steps, file_prefix)

function Molly.log_property!(logger::StateLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        save_reduced_lj_state(s, logger.prefix * "_$(step_n).txt")
    end
end

###reduced temperature logger (kb=1)
struct TemperatureLoggerReduced{T}
    n_steps::Int
    temperatures::Vector{T}
end

TemperatureLoggerReduced(T, n_steps::Integer) = TemperatureLoggerReduced(n_steps, T[])
TemperatureLoggerReduced(n_steps::Integer) = TemperatureLoggerReduced(Float64, n_steps)

function Molly.log_property!(logger::TemperatureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.temperatures, temperature_reduced(s))
    end
end

###virial logger

struct VirialLogger{T}
    n_steps::Int
    energies::Vector{T}
end

VirialLogger(T, n_steps::Integer) = VirialLogger{T}(n_steps, T[])
VirialLogger(n_steps::Integer) = VirialLogger(n_steps, Float64[])


function Molly.log_property!(logger::VirialLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, pair_virial(s, neighbors))
    end
end

###pressure logger
struct PressureLoggerReduced{T}
    n_steps::Int
    pressures::Vector{T}
end

PressureLoggerReduced(T, n_steps::Integer) = PressureLoggerReduced(n_steps, T[])
PressureLoggerReduced(n_steps::Integer) = PressureLoggerReduced(n_steps, Float64[])

function Molly.log_property!(logger::PressureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.pressures, pressure(s, neighbors))
    end
end

struct PressureLoggerNVT{P,K}
    n_steps::Int
    T::K
    pressures::Vector{P}
end

PressureLoggerNVT(T, P, n_steps::Integer) = PressureLoggerNVT(n_steps, T, P[])
PressureLoggerNVT(T, n_steps::Integer) = PressureLoggerNVT(n_steps, T, Float64[])

function Molly.log_property!(logger::PressureLoggerNVT, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        W = pair_virial(s, neighbors)
        V = s.box_size[1] * s.box_size[2] * s.box_size[3]
        push!(logger.pressures, (length(s) * logger.T + W / 3) / V)
    end
end

mutable struct SelfDiffusionLogger{D,T}
    n_steps::Int
    last_coords::T
    self_diffusion_coords::T
    self_diffusion_history::Vector{T}
    record_history::Bool
end

SelfDiffusionLogger(n_steps, initial_coords; record_history::Bool=false) = SelfDiffusionLogger{length(initial_coords[1]),typeof(initial_coords)}(n_steps, deepcopy(initial_coords), zero(initial_coords), typeof(initial_coords)[], record_history)

function Molly.log_property!(logger::SelfDiffusionLogger, s::System, neighbors=nothing, step_n::Integer=0)

    @. logger.self_diffusion_coords += unwrap_coords_vec.(logger.last_coords, s.coords, (s.box_size,)) - logger.last_coords
    @. logger.last_coords = s.coords

    if logger.record_history && (step_n % logger.n_steps == 0)
        push!(logger.self_diffusion_history, deepcopy(logger.self_diffusion_coords))
    end
end

"""
Compute the periodic image of c2 closest to c1 (modulo the side_length), useful for tracking unperiodized coordinates.
"""
function unwrap_coords(c1, c2, side_length)
    if c1 < c2
        return (c1 - c2 + side_length < c2 - c1) ? c2 - side_length : c2
    else
        return (c2 + side_length - c1 < c1 - c2) ? c2 + side_length : c2
    end
end

unwrap_coords_vec(c1, c2, box_size) = unwrap_coords.(c1, c2, box_size)

mutable struct AutoCorrelationLogger{F,T}
    n_steps::Int

    observable::F
    args::T

    corr_length::Int
    Array{}

end


mutable struct CrossCorrelationLogger{F1,F2,T1,T2}
    n_steps

    observable1::F1
    args1::T1

    observable2::F2
    args2::T2

end
