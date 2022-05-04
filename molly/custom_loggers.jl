export 
    HamiltonianLogger,
    KineticEnergyLoggerNoDims,
    StateLogger,
    TemperatureLoggerReduced,
    VirialLogger,
    PressureLoggerReduced,
    PressureLoggerNVT,
    SelfDiffusionLogger,
    GeneralObservableLogger,
    TimeCorrelationLogger,
    TimeCorrelationLoggerVec,
    AutoCorrelationLogger,
    AutoCorrelationLoggerVec,
    ElapsedTimeLogger,
    LogLogger,
    AverageObservableLogger

###Total energy logger

struct HamiltonianLogger{T}
    log_freq::Int
    energies::Vector{T}
end

HamiltonianLogger(T, n_steps::Integer) = HamiltonianLogger(n_steps, T[])
HamiltonianLogger(n_steps::Integer) = HamiltonianLogger(Float32, n_steps)

function Molly.log_property!(logger::HamiltonianLogger, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.energies, Molly.kinetic_energy_noconvert(s) + Molly.potential_energy(s, neighbors))
end

###Dimensionless kinetic energy logger--- Molly version (as of now) requires physical units

struct KineticEnergyLoggerNoDims{T}
    log_freq::Int
    energies::Vector{T}
end

KineticEnergyLoggerNoDims(T, log_freq::Integer) = KineticEnergyLoggerNoDims(log_freq, T[])
KineticEnergyLoggerNoDims(log_freq::Integer) = KineticEnergyLoggerNoDims(Float64, log_freq)

function Molly.log_property!(logger::KineticEnergyLoggerNoDims, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.energies, Molly.kinetic_energy_noconvert(s))
end

###reduced temperature logger (kb=1)
struct TemperatureLoggerReduced{T}
    log_freq::Int
    temperatures::Vector{T}
end

TemperatureLoggerReduced(T, log_freq::Integer) = TemperatureLoggerReduced(log_freq, T[])
TemperatureLoggerReduced(log_freq::Integer) = TemperatureLoggerReduced(Float64, log_freq)

function Molly.log_property!(logger::TemperatureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.temperatures, temperature_reduced(s))
end

###virial logger

struct VirialLogger{T}
    log_freq::Int
    energies::Vector{T}
end

VirialLogger(T, log_freq::Integer) = VirialLogger{T}(log_freq, T[])
VirialLogger(log_freq::Integer) = VirialLogger(log_freq, Float64[])


function Molly.log_property!(logger::VirialLogger, s::System, neighbors=nothing, log_freq::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.energies, pair_virial(s, neighbors))
end

###pressure logger
struct PressureLoggerReduced{T}
    log_freq::Int
    pressures::Vector{T}
end

PressureLoggerReduced(T, log_freq::Integer) = PressureLoggerReduced(log_freq, T[])
PressureLoggerReduced(log_freq::Integer) = PressureLoggerReduced(log_freq, Float64[])

function Molly.log_property!(logger::PressureLoggerReduced, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.pressures, pressure(s, neighbors))
end

struct PressureLoggerNVT{P,K}
    log_freq::Int
    T::K
    pressures::Vector{P}
end

PressureLoggerNVT(T, P, log_freq::Integer) = PressureLoggerNVT(log_freq, T, P[])
PressureLoggerNVT(T, log_freq::Integer) = PressureLoggerNVT(log_freq, T, Float64[])

function Molly.log_property!(logger::PressureLoggerNVT, s::System, neighbors=nothing, step_n::Integer=0)
    if step_n % logger.log_freq == 0
        W = pair_virial(s, neighbors)
        V = s.box_size[1] * s.box_size[2] * s.box_size[3]
        push!(logger.pressures, (length(s) * logger.T + W / 3) / V)
    end
end

mutable struct SelfDiffusionLogger{T}
    last_coords::T
    self_diffusion_coords::T
end

SelfDiffusionLogger(initial_coords) = SelfDiffusionLogger{typeof(initial_coords)}(deepcopy(initial_coords), zero(initial_coords))

function Molly.log_property!(logger::SelfDiffusionLogger, s::System, neighbors=nothing, step_n::Integer=0)
    logger.self_diffusion_coords .+= unwrap_coords_vec.(logger.last_coords, s.coords, (s.box_size,)) - logger.last_coords
    logger.last_coords .= s.coords
end

"""
Compute the periodic image of c2 closest to c1 (modulo the side_length), useful for tracking unperiodized coordinates.
"""
unwrap_coords(c1, c2, side_length)=c2+round((c1-c2)/side_length)*side_length
unwrap_coords_vec(c1, c2, box_size) = unwrap_coords.(c1, c2, box_size)

struct GeneralObservableLogger{T}
    observable::Function
    log_freq::Int64
    history::Vector{T}
end

GeneralObservableLogger(T::DataType, observable::Function, n_steps::Integer) = GeneralObservableLogger{T}(observable, n_steps, T[])
GeneralObservableLogger(observable::Function, n_steps::Integer) = GeneralObservableLogger(Float64, observable, n_steps)

function Molly.log_property!(logger::GeneralObservableLogger, s::System, neighbors=nothing, step_n::Integer=0)
    (step_n % logger.log_freq != 0) && return
    push!(logger.history, logger.observable(s, neighbors))
end


"""
A time correlation logger, which estimates the time correlation function
C(t)=[⟨A(t)B(0)⟩-⟨A(0)B(0)⟩]/√[⟨A(0)²⟩⟨B(0)²⟩]
for observables A and B, which are functions of the form f(sys::System,neighbors=nothing)
the output of A and B can be numeric or vectors (including vectors of SVectors, like the positions or velocities), in which case the products inside the brackets are dot products.
"""
mutable struct TimeCorrelationLogger{T_A,T_sq_A}

    observableA::Function
    observableB::Function

    n_correlation::Integer

    history_A::Vector{T_A}
    history_B::Vector{T_A}

    sum_offset_products::Vector{T_sq_A}

    n_timesteps::Int64

    sum_A::T_A
    sum_B::T_A

    sum_sq_A::T_sq_A
    sum_sq_B::T_sq_A

end


function TimeCorrelationLogger(TA::DataType, observableA::Function, observableB::Function, observable_length::Integer,n_correlation::Integer)
    ini_sum_A= (observable_length>1) ? zeros(TA,observable_length) : zero(TA)
    ini_sum_B= zero(ini_sum_A)
    
    ini_sum_sq_A=dot(ini_sum_A,ini_sum_A)
    ini_sum_sq_B=zero(ini_sum_sq_A)

    T_sum_A=typeof(ini_sum_A)
    T_sum_sq_A=typeof(ini_sum_sq_A)

    return TimeCorrelationLogger{T_sum_A,T_sum_sq_A}(observableA, observableB, n_correlation, T_sum_A[], T_sum_A[], zeros(T_sum_sq_A, n_correlation), 0, ini_sum_A, ini_sum_B, ini_sum_sq_A, ini_sum_sq_B)

end

#Unspecified type defaults to floating-point valued observables
TimeCorrelationLogger(observableA::Function, observableB::Function, n_correlation::Integer) = TimeCorrelationLogger(Float64, observableA, observableB,1,n_correlation)

##Constructors for vector-valued observables
TimeCorrelationLoggerVec(N_atoms::Integer, dim::Integer, T::DataType, observableA::Function, observableB::Function, n_correlation::Integer) = (dim>1) ? TimeCorrelationLogger(SVector{dim,T}, observableA, observableB,N_atoms, n_correlation) : TimeCorrelationLogger(T, observableA, observableB,N_atoms, n_correlation)
TimeCorrelationLoggerVec(N_atoms::Integer, dim::Integer, observableA::Function, observableB::Function, n_correlation::Integer) = TimeCorrelationLoggerVec(N_atoms, dim, Float64, observableA, observableB, n_correlation)

## Convenience constructors for autocorrelations
AutoCorrelationLogger(T::DataType, observable::Function, n_correlation::Integer) = TimeCorrelationLogger(T, observable, observable,observable_length, 1,n_correlation)
AutoCorrelationLogger(observable::Function, n_correlation::Integer) = TimeCorrelationLogger(Float64, observable, n_correlation)
AutoCorrelationLoggerVec(N_atoms::Integer, dim::Integer, T::DataType, observable::Function, n_correlation::Integer) = TimeCorrelationLoggerVec(N_atoms, dim, T,observable, observable, n_correlation)
AutoCorrelationLoggerVec(N_atoms::Integer, dim::Integer, observable::Function, n_correlation::Integer) = AutoCorrelationLoggerVec(N_atoms, dim, Float64, observable, n_correlation)


function Molly.log_property!(logger::TimeCorrelationLogger, s::System, neighbors=nothing, step_n::Integer=0)

    #compute observables
    A = logger.observableA(s, neighbors)
    B = (logger.observableA!=logger.observableB) ? logger.observableB(s, neighbors) : A

    logger.n_timesteps += 1

    #update history lists
    (logger.n_timesteps > logger.n_correlation) && (popfirst!(logger.history_A); popfirst!(logger.history_B))

    push!(logger.history_A, A)
    push!(logger.history_B, B)

    #update running averages (numerically stable method)
    logger.sum_A += A
    logger.sum_B += B

    logger.sum_sq_A += dot(A, A)
    logger.sum_sq_B += dot(B, B)


    B1 = first(logger.history_B)

    for (i, Ai) = enumerate(logger.history_A)
        n_samples = logger.n_timesteps - i + 1
        if n_samples > 0
            logger.sum_offset_products[i] += dot(Ai, B1)
        end
    end

end



struct ElapsedTimeLogger
    initial_time::Dates.DateTime
end

ElapsedTimeLogger()=ElapsedTimeLogger(Dates.now())
Molly.log_property!(logger::ElapsedTimeLogger,S::System,neighbors=nothing,step_n::Integer=0)=nothing

struct LogLogger
    logger_table::Vector{Tuple{Symbol,String,Int64,Bool,String}} 
end

function LogLogger(logger_names::Vector{Symbol},logging_files::Vector{String},logging_freqs::Vector{Int64},decorate_outputs::Vector{Bool},write_modes::Vector{String})
    return LogLogger([zip(logger_names,logging_files,logging_freqs,decorate_outputs,write_modes)...])
end

function Molly.log_property!(logger::LogLogger,s::System,neighbors=nothing,step_n::Integer=0)
    for (name,file,freq,decorate,mode)=logger.logger_table
        if step_n%freq==0
            f=open(file,mode)
            (decorate) && println(f,"---- Step : $(step_n) ---- Logger : $(name) ----")
            log_to_file!(s.loggers[name],f)
            close(f)
        end
    end
end

function log_to_file!(logger::GeneralObservableLogger,file::IOStream)
    write(file,logger.history)
    empty!(logger.history)
end

function log_to_file!(logger::TimeCorrelationLogger,file::IOStream)#write input utility
    write(file,length(logger.sum_A))
    write(file,logger.sum_A)
    write(file,logger.sum_B)
    write(file,logger.sum_sq_A)
    write(file,logger.sum_sq_B)
    write(file,length(logger.sum_offset_products))
    write(file,logger.sum_offset_products)
    write(file,logger.n_timesteps)
end

function log_to_file!(logger::ElapsedTimeLogger,file::IOStream)
    println(file,Dates.canonicalize(Dates.now()-logger.initial_time))
end

function log_to_file!(logger::SelfDiffusionLogger,file::IOStream)
    write(file,logger.self_diffusion_coords)
end

###State logger (writes state of system to external file --- NOT GENERAL)

struct StateLogger
    s::System
end

Molly.log_property!(logger::StateLogger,s::System,neighbors=nothing,step_n::Integer=0)=nothing

function log_to_file!(logger::StateLogger, file::IOStream)
    save_reduced_state(logger.s, file)
end

mutable struct AverageObservableLogger{T}
    observable::Function
    log_freq::Int64
    sum::T
    n_samples::Int64
end

AverageObservableLogger(T::DataType,observable::Function,log_freq::Integer=1)=AverageObservableLogger{T}(observable,log_freq,zero(T),0)
AverageObservableLogger(observable::Function,log_freq::Integer=1)=AverageObservableLogger(Float64,observable,log_freq)

function Molly.log_property!(logger::AverageObservableLogger,s::System,neighbors=nothing,step_n::Integer=0)

    if step_n % logger.log_freq ==0
        O=logger.observable(s,neighbors)
        logger.n_samples+=1
        logger.sum+=O
    end

end

function log_to_file!(logger::AverageObservableLogger,file::IOStream)
    write(file,logger.n_samples)
    write(file,logger.sum/logger.n_samples)
end