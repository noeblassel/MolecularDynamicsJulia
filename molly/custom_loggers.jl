###Total energy logger

struct HamiltonianLogger{T}
    n_steps::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer) = HamiltonianLogger(n_steps, T[])
HamiltonianLogger(n_steps::Integer) = HamiltonianLogger(Float32, n_steps)

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors = nothing, step_n::Integer = 0)
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

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors = nothing, step_n::Integer = 0)
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

function Molly.log_property!(logger::StateLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        save_reduced_lj_state(s, logger.prefix * "_$(step_n).txt")
    end
end

###reduced temperature logger (kb=1)
struct TemperatureLoggerReduced{T}
    n_steps::Int
    temperatures::Vector{T}
end

TemperatureLoggerReduced(T, n_steps::Integer)=TemperatureLoggerReduced(n_steps,T[])
TemperatureLoggerReduced(n_steps::Integer)=TemperatureLoggerReduced(Float64, n_steps)

function Molly.log_property!(logger::TemperatureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.temperatures,temperature_reduced(s))
    end
end

###virial logger

struct VirialLogger{T}
    n_steps::Int
    energies::Vector{T}
end

VirialLogger(T, n_steps::Integer) = VirialLogger{T}(n_steps, T[])
VirialLogger(n_steps::Integer) = VirialLogger(n_steps, Float64[])


function Molly.log_property!(logger::VirialLogger, s::System, neighbors = nothing, step_n::Integer = 0)
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

function Molly.log_property!(logger::PressureLoggerReduced, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        push!(logger.pressures, pressure(s,neighbors))
    end
end