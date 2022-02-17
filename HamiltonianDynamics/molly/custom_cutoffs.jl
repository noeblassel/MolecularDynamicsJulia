"""The spline cutoff, yielding a continuously differentiable intepolation between the true potential at radius activation_dist and zero at distance dist_cutoff
"""
struct CubicSplineCutoff{D,S,I}

    dist_cutoff::D
    sqdist_cutoff::S
    inv_sqdist_cutoff::I

    activation_dist::S#Actually the square of the activation distance
    inv_activation_dist::I#Actually the inverse square activation distance
    activation_dist_true::D

    delta::D

end

function CubicSplineCutoff(activation_dist, dist_cutoff)
    dc = dist_cutoff
    ad = activation_dist^2 ##The problem arises from variable naming at lennard_jones.jl l. 71 (function force(::LennardJones,args...))
    δ = dc - activation_dist

    if δ <= 0
        error("The cutoff radius must be strictly larger than the activation radius.")
    end

    (D, S, I) = typeof.([dc, dc^2, inv(dc^2)])

    return CubicSplineCutoff{D,S,I}(dc, dc^2, inv(dc^2), ad, inv(ad), sqrt(ad), δ)
end

Molly.cutoff_points(::Type{CubicSplineCutoff{D,S,I}}) where {D,S,I} = 2


@fastmath function Molly.force_divr_cutoff(cutoff::CubicSplineCutoff, r2, inter, params)
    r = √r2
    rs = cutoff.activation_dist_true
    t = (r - rs) / cutoff.delta

    ps = Molly.potential(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params)
    fs = -Molly.force_divr_nocutoff(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params) * rs
    ##cubic spline interpolation derivative
    return (6t^2 - 6t) * ps / cutoff.delta + (3t^2 - 4t + 1) * fs

end

@fastmath function Molly.potential_cutoff(cutoff::CubicSplineCutoff, r2, inter, params)
    r = √r2
    rs = cutoff.activation_dist_true
    t = (r - rs) / cutoff.delta

    ps = Molly.potential(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params)
    fs = -Molly.force_divr_nocutoff(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params) * rs
    ##cubic spline interpolation
    return (2t^3 - 3t^2 + 1) * ps + (t^3 - 2t^2 + t) * cutoff.delta * fs
end


#Shifted force cutoff is broken

"""
    ShiftedForceCutoff_(dist_cutoff)

Cutoff that shifts the force to be continuous at a specified cutoff point.
"""
struct ShiftedForceCutoff_{D, S, I}
    dist_cutoff::D
    sqdist_cutoff::S
    inv_sqdist_cutoff::I
end

function Molly.ShiftedForceCutoff_(dist_cutoff)
    return ShiftedForceCutoff_(dist_cutoff, dist_cutoff ^ 2, inv(dist_cutoff ^ 2))
end

Molly.cutoff_points(::Type{ShiftedForceCutoff_{D, S, I}}) where {D, S, I} = 1

function Molly.force_divr_cutoff(cutoff::ShiftedForceCutoff_, r2, inter, params)
    return Molly.force_divr_nocutoff(inter, r2, inv(r2), params) - Molly.force_divr_nocutoff(
                                inter, cutoff.sqdist_cutoff, cutoff.inv_sqdist_cutoff, params)
end

@fastmath function Molly.potential_cutoff(cutoff::ShiftedForceCutoff_, r2, inter, params)
    invr2 = inv(r2)
    r = √r2
    rc = cutoff.dist_cutoff
    fc = Molly.force_divr_nocutoff(inter, cutoff.sqdist_cutoff, cutoff.inv_sqdist_cutoff, params) * r

    Moly.potential(inter, r2, invr2, params) - (r - rc) * fc -
        Molly.potential(inter, cutoff.sqdist_cutoff, cutoff.inv_sqdist_cutoff, params)
end