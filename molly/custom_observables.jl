export
    pair_virial,
    pressure,
    temperature_reduced,
    long_range_virial_correction,
    MobilityObservable

function pair_virial(s::System, neighbors = nothing, lrc = false)

    W = 0.0 * s.energy_units
    N=length(s)

    for inter = s.pairwise_inters
        if isnothing(neighbors) || (!inter.nl_only) #if neighbor list is not given or unused, then double loop
            for i = 1:N
                for j = 1:i-1

                    σ = inter.lorentz_mixing ? (s.atoms[i].σ + s.atoms[j].σ) / 2 : sqrt(s.atoms[i].σ * s.atoms[j].σ)
                    ϵ = sqrt(s.atoms[i].ϵ * s.atoms[j].ϵ)

                    r2 = Molly.square_distance(i, j, s.coords, s.box_size)

                    cutoff = inter.cutoff
                    C = typeof(cutoff)
                    σ2 = σ^2
                    params = (σ2, ϵ)

                    f_divr = 0.0

                    #this piece of logic could be packaged into a function

                    if Molly.cutoff_points(C) == 0
                        f_divr = Molly.force_divr_nocutoff(inter, r2, inv(r2), params)
                    elseif Molly.cutoff_points(C) == 1
                        f_divr = r2 > cutoff.sqdist_cutoff ? 0.0 * inter.force_units : Molly.force_divr_cutoff(cutoff, r2, inter, params)
                    elseif Molly.cutoff_points(C) == 2
                        if r2 > cutoff.sqdist_cutoff
                            f_divr = 0.0 * inter.force_units

                        elseif r2 < cutoff.sqdist_activation
                            f_divr = Molly.force_divr_nocutoff(inter, r2, inv(r2), params)
                        else
                            f_divr = Molly.force_divr_cutoff(cutoff, r2, inter, params)
                        end
                    end

                    W += f_divr * r2
                end
            end


        else
            for ni = 1:neighbors.n
                i, j, _ = neighbors.list[ni]
                σ = inter.lorentz_mixing ? (s.atoms[i].σ + s.atoms[j].σ) / 2 : sqrt(s.atoms[i].σ * s.atoms[j].σ)
                ϵ = sqrt(s.atoms[i].ϵ * s.atoms[j].ϵ)

                r2 = Molly.square_distance(i, j, s.coords, s.box_size)

                cutoff = inter.cutoff
                C = typeof(cutoff)
                σ2 = σ^2
                params = (σ2, ϵ)

                f_divr = 0.0

                #this piece of logic could be packaged into a function

                if Molly.cutoff_points(C) == 0
                    f_divr = Molly.force_divr_nocutoff(inter, r2, inv(r2), params)
                elseif Molly.cutoff_points(C) == 1
                    f_divr = r2 > cutoff.sqdist_cutoff ? 0.0 * inter.force_units : Molly.force_divr_cutoff(cutoff, r2, inter, params)
                elseif Molly.cutoff_points(C) == 2
                    if r2 > cutoff.sqdist_cutoff
                        f_divr = 0.0 * inter.force_units

                    elseif r2 < cutoff.sqdist_activation
                        f_divr = Molly.force_divr_nocutoff(inter, r2, inv(r2), params)
                    else
                        f_divr = Molly.force_divr_cutoff(cutoff, r2, inter, params)
                    end
                end

                W += f_divr * r2
            end
        end

        if lrc
            W += long_range_virial_correction(s, inter)
        end
    end
    return W
end


#not unit safe

function pressure(s::System, neighbors = nothing)
    l1, l2, l3 = s.box_size
    V = l1 * l2 * l3
    K = Molly.kinetic_energy_noconvert(s)
    W = pair_virial(s, neighbors)
    
    return (2K + W) / 3V
end

function temperature_reduced(s::System)##ie when kb=1
    N = length(s)
    ke = Molly.kinetic_energy_noconvert(s)

    return 2ke / (3N - 3)
end

function long_range_virial_correction(s::System, inter::LennardJones)
        !hasproperty(inter, :cutoff) && return 0.0 * s.energy_units
    isa(inter.cutoff, NoCutoff) && return 0.0 * s.energy_units


    l1, l2, l3 = s.box_size
    V = l1 * l2 * l3
    N = length(s)

    r_cm3 = inv(inter.cutoff.dist_cutoff^3)
    
end


function MobilityObservable(ff::Vector)
    R(s::System,neighbors=nothing)=sum(dot.(ff,s.velocities))
    return R
end