export
    NonGradientForceNEMD,
    OneDriftNEMD,
    TwoDriftNEMD,
    ColorDriftNEMD

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
    return NonGradientForceNEMD{D,T}(N,η,ff)
end

function TwoDriftNEMD(N::Integer,η::Real;D::Integer=3,T::DataType=Float64)
    ff=zeros(SVector{D,T},N)
    ff[1]=SVector{D,T}(vcat(one(T),zeros(T,D-1)))
    ff[2]=SVector{D,T}(vcat(-one(T),zeros(T,D-1)))
    normalize!(ff)
    return NonGradientForceNEMD{D,T}(N,η,ff)
end

function ColorDriftNEMD(N::Integer,η::Real;D::Integer=3,T::DataType=Float64)
    ff=[SVector{D,T}(vcat(one(T)*(-1.0) ^(i+1),zeros(T,D-1))) for i=1:N]
    normalize!(ff)
    return NonGradientForceNEMD{D,T}(N,η,ff)
end