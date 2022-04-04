export 
    place_atoms_on_3D_lattice,
    init_velocities,
    reduced_velocity


function place_atoms_on_3D_lattice(N_per_dim::Integer,box_size)
    (Lx,Ly,Lz)=box_size
    reshape([SVector(i * Lx/N_per_dim, j * Ly/N_per_dim, k * Lz/N_per_dim) for i = 1:N_per_dim, j = 1:N_per_dim, k = 1:N_per_dim],N_per_dim^3)
end

function init_velocities(T::Ftype,M::Vector{Ftype},k::Ftype) where {Ftype<:Real}
    N=length(M)
    return [reduced_velocity(T,M[i],k) for i=1:N]
end

reduced_velocity_lj(T::Real,M::Real,k::Real) = sqrt(k*T/M) * (@SVector randn(3))
reduced_velocity(T::Real,M::Real,k::Real) = sqrt(k*T/M) * (SVector{3,Float64}(randn(3)))