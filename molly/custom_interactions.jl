export
    NonGradientForceNEMD,
    OneDriftNEMD,
    TwoDriftNEMD,
    ColorDriftNEMD,
    TransverseForceProfile,
    PiecewiseConstantForceProfile,
    PiecewiseLinearForceProfile,
    SinusoidalForceProfile

struct NonGradientForceNEMD{D,T}
    N::Integer
    η::Real
    force_field::Vector{SVector{D,T}}
end

Molly.forces(inter::NonGradientForceNEMD,s::System,neighbors)=inter.η*inter.force_field

## specific forcings, always such that |F|=1
function OneDriftNEMD(N::Integer,η::Real;D::Integer=3,T::DataType=Float64)
    ff=zeros(SVector{D,T},N)
    ff[1]=SVector{D,T}(vcat(one(T),zeros(T,D-1)))
    normalize!(ff)
    #shuffle!(ff)
    return NonGradientForceNEMD{D,T}(N,η,ff)
end

function TwoDriftNEMD(N::Integer,η::Real;D::Integer=3,T::DataType=Float64)
    ff=zeros(SVector{D,T},N)
    ff[1]=SVector{D,T}(vcat(one(T),zeros(T,D-1)))
    ff[2]=SVector{D,T}(vcat(-one(T),zeros(T,D-1)))
    normalize!(ff)
    shuffle!(ff)
    return NonGradientForceNEMD{D,T}(N,η,ff)
end

function ColorDriftNEMD(N::Integer,η::Real;D::Integer=3,T::DataType=Float64)
    ff=[SVector{D,T}(vcat(one(T)*(-1.0) ^(i+1),zeros(T,D-1))) for i=1:N]
    normalize!(ff)
    #shuffle!(ff)
    return NonGradientForceNEMD{D,T}(N,η,ff)
end

struct TransverseForceProfile{TF}
    ξ::Real
    F_no_units::TF # a function which maps the a single position coordinate to a (unitless) force
end

TransverseForceProfile(ξ::Real,F::Function)=TransverseForceProfile{typeof(F)}(ξ,F)
x_coord_to_vec(val::T,dim::Integer) where {T} = SVector{dim,T}(vcat(val,zeros(T,dim-1)))

function Molly.forces(inter::TransverseForceProfile,s::System{D},neighbors) where {D}
    y_coords= getindex.(s.coords,(2,))
    force_vec=inter.F_no_units.(y_coords)
    return inter.ξ * x_coord_to_vec.(force_vec,(D,)) * s.force_units #units are handled here
end

function SinusoidalForceProfile(;ξ::Real,L::T) where {T}
    sinusoidal_force(y::T) = sin(2π*ustrip(y/L))
    return TransverseForceProfile(ξ,sinusoidal_force)
end

function PiecewiseConstantForceProfile(;ξ::Real,L::T) where {T}
    piecewise_constant_force(y::T) = ( y < L / 2 ) ? -1.0 : 1.0
    return TransverseForceProfile(ξ,piecewise_constant_force)
end

function PiecewiseLinearForceProfile(;ξ::Real,L::T) where {T}
    piecewise_linear_force(y::T) = ( y < L / 2 ) ? ustrip(4*(y-L/4)/L) : ustrip(4*(3L/4-y)/L)
    return TransverseForceProfile(ξ,piecewise_linear_force)
end
