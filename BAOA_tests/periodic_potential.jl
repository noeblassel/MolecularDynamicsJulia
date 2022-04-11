struct SinusoidalPotential <: Interaction
end

function Molly.potential_energy(inter::SinusoidalPotential,s::System,neighbors=nothing)
    V=zero(first(s.coords))

    for q in s.coords
        @. V+=sin(2π*q/s.box_size)
    end

    return sum(V)*s.energy_units
end

function Molly.forces(inter::SinusoidalPotential,s::System,neighbors=nothing)
    fs=zero(s.coords)

    for i=1:length(s)
    fs[i]-=(2π./s.box_size).*cos.(2π*s.coords[i]./s.box_size)
    end
    return fs*s.force_units
end
