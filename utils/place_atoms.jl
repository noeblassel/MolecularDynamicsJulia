function place_atoms_on_lattice(N_per_dim::Integer,box_size)
    (lx,ly,lz)=box_size/N_per_dim
    ixs=[(i,j,k) for i=1:N_per_dim,j=1:N_per_dim,k=1:N_per_dim]

    return MVector{N_per_dim^3}(SVector(i*lx,j*ly,k*lz) for (i,j,k)=ixs)
end

function init_velocities(T::Ftype,M::Vector{Ftype},k::Ftype) where {Ftype<:Real}
    N=length(M)
    return MVector{N}(reduced_velocity(T,M[i],k) for i=1:N)
end

reduced_velocity(T::Real,M::Real,k::Real) = sqrt(k*T/M) * (@SVector randn(3))