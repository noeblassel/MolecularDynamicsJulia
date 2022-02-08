struct PressureLoggerLJ{T}
    n_steps::Int
    pressures::Vector{T}
end

PressureLoggerLJ(T, n_steps::Integer)=PressureLoggerLJ(n_steps,T[])
PressureLoggerLJ(n_steps::Integer)=PressureLoggerLJ(Float32, n_steps)

function Base.show(io::IO, pl::PressureLoggerLJ)
    print(io, "PressureLoggerLJ{", eltype(eltype(pl.pressures)), "} with n_steps ",
            pl.n_steps, ", ", length(pl.pressures), " pressures recorded")
end

function Molly.log_property!(logger::PressureLoggerLJ, s::System, neighbors=nothing, step_n::Integer=0)
    N=length(s.coords)
    if step_n % logger.n_steps == 0
        ke=Molly.kinetic_energy_noconvert(s)
        W=0
        for i=1:N
            for j=1:i-1
                r_ij2=Molly.square_distance(i,j,s.coords,s.box_size)
                W-=Molly.force_divr_nocutoff(s.general_inters[1],r_ij2,inv(r_ij2),(1.0,1.0))*r_ij2 #force_divr_nocutoff(LJ,r2,1/r2,(σ,ϵ))=-LJ'(r)/r
            end
        end

        P=2*ke+W
        P/=3*s.box_size[1]^3
        push!(logger.pressures,P)
    end
end