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

    Vs = Molly.potential(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params)
    dVs = -Molly.force_divr_nocutoff(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params) * rs
    ##cubic spline interpolation derivative
    return -((6t^2 - 6t) * Vs / cutoff.delta + (3t^2 - 4t + 1) * dVs)/r

end

@fastmath function Molly.potential_cutoff(cutoff::CubicSplineCutoff, r2, inter, params)
    r = √r2
    rs = cutoff.activation_dist_true
    t = (r - rs) / cutoff.delta

    Vps = Molly.potential(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params)
    dVs = -Molly.force_divr_nocutoff(inter, cutoff.activation_dist, cutoff.inv_activation_dist, params) * rs
    ##cubic spline interpolation
    return (2t^3 - 3t^2 + 1) * Vs + (t^3 - 2t^2 + t) * cutoff.delta * dVs
end