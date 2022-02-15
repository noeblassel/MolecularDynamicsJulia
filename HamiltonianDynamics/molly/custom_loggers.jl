###Total energy logger

struct HamiltonianLogger{T}
    n_steps::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer)=HamiltonianLogger(n_steps,T[])
HamiltonianLogger(n_steps::Integer)=HamiltonianLogger(Float32, n_steps)

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies,Molly.kinetic_energy_noconvert(s)+Molly.potential_energy(s,neighbors))
    end
end

###Dimensionless kinetic energy logger--- Molly version requires physical units

struct KineticEnergyLoggerNoDims{T}
    n_steps::Int
    energies::Vector{T}
end

KineticEnergyLoggerNoDims(T, n_steps::Integer)=KineticEnergyLoggerNoDims(n_steps,T[])
KineticEnergyLoggerNoDims(n_steps::Integer)=KineticEnergyLoggerNoDims(Float64, n_steps)

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        push!(logger.energies,Molly.kinetic_energy_noconvert(s))
    end
end

###State logger (writes state of system to external file)

struct StateLogger
    n_steps::Int
    prefix::AbstractString
end

StateLogger(n_steps::Integer) = StateLogger(n_steps, "logfile")
StateLogger(n_steps::Integer, file_prefix::AbstractString) = StateLogger(n_steps, file_prefix)

function Molly.log_property!(logger::StateLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        save_reduced_lj_state(s,logger.prefix*"_$(step_n).txt")
    end
end

###redundant (TemperatureLogger is fixed-- GH version is not up to date )

"""###reduced temperature logger (kb=1)
struct ReducedTemperatureLogger{T}
    n_steps::Int
    temperatures::Vector{T}
end

ReducedTemperatureLogger(T, n_steps::Integer)=ReducedTemperatureLogger(n_steps,T[])
ReducedTemperatureLogger(n_steps::Integer)=ReducedTemperatureLogger(Float64, n_steps)

function Molly.log_property!(logger::ReducedTemperatureLogger, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.n_steps == 0
        N=length(s)
        ke=Molly.kinetic_energy_noconvert(s)
        T=2ke/(3N)
        push!(logger.temperatures,T)
    end
end"""

###virial logger

struct VirialLogger{T}
    n_steps::Int
    energies::Vector{T}
end

VirialLogger(T, n_steps::Integer) = VirialLogger{T}(n_steps, T[])
VirialLogger(n_steps::Integer) = VirialLogger(n_steps, Float64[])


function Molly.log_property!(logger::VirialLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, virial(s,neighbors))
    end
end


function virial(s::System,neighbors=nothing)
    r_c2=s.general_inters[1].cutoff.sqdist_cutoff
    N = length(s)
    W = 0
    for i = 1:N
        for j = 1:i-1
            r_ij2 = Molly.square_distance(i, j, s.coords, s.box_size)
            if r_ij2 < r_c2
                W += Molly.force_divr_cutoff(s.general_inters[1].cutoff, r_ij2, s.general_inters[1], (1.0, 1.0)) * r_ij2 #force_divr_nocutoff(LJ,r2,1/r2,(σ,ϵ))=-LJ'(r)/r
            end
        end
    end
    return W
end

###pressure logger
struct PressureLoggerReduced{T}
    n_steps::Int
    pressures::Vector{T}
end

PressureLoggerReduced(T, n_steps::Integer) = PressureLoggerReduced(n_steps, T[])
PressureLoggerReduced(n_steps::Integer) = PressureLoggerReduced(n_steps, Float64[])
PressureLoggerReduced(n_steps::Integer) = PressureLoggerReduced(n_steps, Float64[])


function Molly.log_property!(logger::PressureLoggerReduced, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        N = length(s)
        l1,l2,l3=s.box_size
        V=l1*l2*l3
        K = Molly.kinetic_energy_noconvert(s)
        W=virial(s,neighbors)
        P=(2K+W)/3V
        push!(logger.pressures, P)
    end
end