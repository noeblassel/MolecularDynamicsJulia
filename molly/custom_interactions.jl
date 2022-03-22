struct SingleDriftNEMD{D,T} 
    N_atoms::Integer
    ix::Integer

    η::Real
    F::SVector{D,T}

    force_field::Vector{SVector{D,T}}
end

function SingleDriftNEMD(N_atoms::Integer,ix::Integer,η::Real,F::SVector{D,T}) where {D,T}
    ff=zeros(T=typeof(F),dims=N_atoms)
    ff[ix]=η*F
    SingleDriftNEMD{D,T}(N_atoms,ix,η,F ,ff)

end

SingleDriftNEMD(η::Float64,ix::Int64,F::Float64)=SingleDriftNEMD{3,Float64}(η,ix,SVector(F,F,F))
SingleDriftNEMD(η::Float64,F)=SingleDriftNEMD(η,0,F)

forces(inter, s::System{D}, neighbors) where {D}
    return inter.force_field ./ mass.(s.atoms)
end


struct ColorDriftNEMD{D,T}
    N_atoms::Int64
    ix::Int64

    η::Real
    F::SVector{D,T}

    force_field::Vector{SVector{D,T}}
end

ColorDriftNEMD(η::Float64,F::Float64)=ColorDriftNEMD{3,Float64}(η,SVector(F,F,F))
ColorDriftNEMD(η::Float64,F::SVector)=ColorDriftNEMD{length(F),typeof(F[1])}(η,SVector(F,F,F))

forces(inter, s::System{D}, neighbors) where {D}
    return inter.force_field ./ mass.(s.atoms)
end