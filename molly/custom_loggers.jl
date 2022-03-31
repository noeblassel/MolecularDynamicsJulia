###Total energy logger

struct HamiltonianLogger{T}
    n_steps::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer) = HamiltonianLogger(n_steps, T[])
HamiltonianLogger(n_steps::Integer) = HamiltonianLogger(Float32, n_steps)

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.energies, Molly.kinetic_energy_noconvert(s) + Molly.potential_energy(s, neighbors))
end

###Dimensionless kinetic energy logger--- Molly version (as of now) requires physical units

struct KineticEnergyLoggerNoDims{T}
    n_steps::Int
    energies::Vector{T}
end

KineticEnergyLoggerNoDims(T, n_steps::Integer) = KineticEnergyLoggerNoDims(n_steps, T[])
KineticEnergyLoggerNoDims(n_steps::Integer) = KineticEnergyLoggerNoDims(Float64, n_steps)

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.energies, Molly.kinetic_energy_noconvert(s))
end

###State logger (writes state of system to external file --- NOT GENERAL)

struct StateLogger
    n_steps::Int
    prefix::AbstractString
end

StateLogger(n_steps::Integer) = StateLogger(n_steps, "logfile")
StateLogger(n_steps::Integer, file_prefix::AbstractString) = StateLogger(n_steps, file_prefix)

function Molly.log_property!(logger::StateLogger, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    save_reduced_lj_state(s, logger.prefix * "_$(step_n).txt")
end

###reduced temperature logger (kb=1)
struct TemperatureLoggerReduced{T}
    n_steps::Int
    temperatures::Vector{T}
end

TemperatureLoggerReduced(T, n_steps::Integer) = TemperatureLoggerReduced(n_steps, T[])
TemperatureLoggerReduced(n_steps::Integer) = TemperatureLoggerReduced(Float64, n_steps)

function Molly.log_property!(logger::TemperatureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.temperatures, temperature_reduced(s))
end

###virial logger

struct VirialLogger{T}
    n_steps::Int
    energies::Vector{T}
end

VirialLogger(T, n_steps::Integer) = VirialLogger{T}(n_steps, T[])
VirialLogger(n_steps::Integer) = VirialLogger(n_steps, Float64[])


function Molly.log_property!(logger::VirialLogger, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.energies, pair_virial(s, neighbors))
end

###pressure logger
struct PressureLoggerReduced{T}
    n_steps::Int
    pressures::Vector{T}
end

PressureLoggerReduced(T, n_steps::Integer) = PressureLoggerReduced(n_steps, T[])
PressureLoggerReduced(n_steps::Integer) = PressureLoggerReduced(n_steps, Float64[])

function Molly.log_property!(logger::PressureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.pressures, pressure(s, neighbors))
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

mutable struct SelfDiffusionLogger{T}
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

    (step_n % logger.n_steps != 0) && return
    push!(logger.self_diffusion_history, deepcopy(logger.self_diffusion_coords))
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

struct GeneralObservableLogger{T}
    observable::Function
    n_steps::Int64
    history::Vector{T}
end

GeneralObservableLogger(T::DataType,observable::Function,n_steps::Integer)=GeneralObservableLogger{T}(observable,n_steps,T[])
GeneralObservableLogger(observable::Function,n_steps::Integer)=GeneralObservableLogger(Float64,observable,n_steps)

function Molly.log_property!(logger::GeneralObservableLogger,s::System,neighbors=nothing,step_n::Integer=0)
    (step_n % logger.n_steps != 0) && return
    push!(logger.history,logger.observable(s,neighbors))
end

mutable struct TimeCorrelationLogger{T}

    observableA::Function
    observableB::Function

    n_correlation::Integer
    
    history_A::Vector{T}
    history_B::Vector{T}

    uncentered_correlations::Vector{T}
    correlations::Vector{T}

    n_timesteps::Int64

    avg_A::T
    avg_B::T
end

TimeCorrelationLogger(T::DataType,observableA::Function,observableB::Function,n_correlation::Integer)=TimeCorrelationLogger{T}(observableA,observableB,n_correlation,T[],T[],zeros(T,n_correlation),zeros(T,n_correlation),0,zero(T),zero(T))
TimeCorrelationLogger(observableA::Function,observableB::Function,n_correlation::Integer)=TimeCorrelationLogger(Float64,observableA,observableB,n_correlation)

AutoCorrelationLogger(T::DataType,observable::Function,n_correlation::Integer)=TimeCorrelationLogger(T,observable,observable,n_correlation)
AutoCorrelationLogger(observable::Function,n_correlation::Integer)=TimeCorrelationLogger(Float64,observable,observable,n_correlation)


function Molly.log_property!(logger::TimeCorrelationLogger,s::System,neighbors=nothing,step_n::Integer=0)
    A=logger.observableA(s,neighbors)
    B=logger.observableB(s,neighbors)

    logger.n_timesteps+=1

    logger.avg_A=((logger.n_timesteps-1)*logger.avg_A+A)/logger.n_timesteps
    logger.avg_B=((logger.n_timesteps-1)*logger.avg_B+B)/logger.n_timesteps

    push!(logger.history_A,A)
    push!(logger.history_B,B)

    (length(logger.history_A) > logger.n_correlation) && (popfirst!(logger.history_A) ; popfirst!(logger.history_B))
    B1=first(logger.history_B)

    for (i,Ai)=enumerate(logger.history_A)
        n_sampled=logger.n_timesteps-i

        if n_sampled>=0
            logger.uncentered_correlations[i]=(logger.uncentered_correlations[i]*n_sampled + Ai*B1)/(n_sampled+1)
        end
    end

    @. logger.correlations=logger.uncentered_correlations - logger.avg_A*logger.avg_B    
end