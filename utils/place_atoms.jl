function place_atoms_on_lattice(N_per_dim::Integer,box_size)
    (Lx,Ly,Lz)=box_size
    reshape([SVector(i * Lx/N_per_dim, j * Ly/N_per_dim, k * Lz/N_per_dim) for i = 1:N_per_dim, j = 1:N_per_dim, k = 1:N_per_dim],N_per_dim^3)
end

reduced_velocity_lj(T::Real) = sqrt(T) * (@SVector randn(3))