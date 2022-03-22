struct SingleDriftNEMD{D,T} 
    N_atoms::Integer
    ix::Integer

    η::Real
    F::SVector{D,T}

    force_field::Vector{SVector{D,T}}
end

function SingleDriftNEMD(N_atoms::Integer,ix::Integer,η::Real,F::SVector{D,T}) where {D,T}
    ff=zeros(SVector{D,T},N_atoms)
    ff[ix]=η*F
    return SingleDriftNEMD{D,T}(N_atoms,ix,η,F ,ff)

end

SingleDriftNEMD(N_atoms::Integer,η::Float64,ix::Integer) where {T}=SingleDriftNEMD(N_atoms,0,η,SVector(1.0,0.0,0.0))

forces(inter::SingleDriftNEMD, s::System{D}, neighbors)=inter.force_field


struct ColorDriftNEMD{D,T}
    N_atoms::Int64

    η::Real
    F::SVector{D,T}

    force_field::Vector{SVector{D,T}}
end

function ColorDriftNEMD(N_atoms::Integer,η::Real,F::SVector{D,T}) where {D,T}
    ff=zeros(SVector{D,T},N_atoms)

    for ix=1:2:length(ff)
        ff[ix]=η*F
    end

    for ix=2:2:length(ff)
        ff[ix]=-η*F
    end

    return ColorDriftNEMD{D,T}(N_atoms,η,F,ff)

end
ColorDriftNEMD(N_atoms::Integer,η::Float64)=ColorDriftNEMD{3,Float64}(η,SVector(1.0,0.0,0.0))

forces(inter::ColorDriftNEMD, s::System, neighbors)=inter.force_field