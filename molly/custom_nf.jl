export CellListMapNeighborFinder2D

mutable struct CellListMapNeighborFinder2D{N, T} <: AbstractNeighborFinder
    nb_matrix::BitArray{2}
    matrix_14::BitArray{2}
    n_steps::Int
    dist_cutoff::T
    # Auxiliary arrays for multi-threaded in-place updating of the lists
    cl::CellListMap.CellList{N, T}
    aux::CellListMap.AuxThreaded{N, T}
    neighbors_threaded::Vector{NeighborList}
end

function CellListMapNeighborFinder2D(;
    nb_matrix,
    matrix_14=falses(size(nb_matrix)),
    n_steps=10,
    x0=nothing,
    unit_cell=nothing,
    number_of_batches=(0, 0), # (0, 0): use default heuristic
    dist_cutoff::T) where T
np = size(nb_matrix, 1)
if isnothing(unit_cell)
if unit(dist_cutoff) == NoUnits
side = max(2 * dist_cutoff, (np * 0.01) ^ (1 / 2))
else
side = max(2 * dist_cutoff, uconvert(unit(dist_cutoff), (np * 0.01u"nm^3") ^ (1 / 2)))
end
sides = SVector(side, side, side)
box = CellListMap.Box(sides, dist_cutoff)
else
box = CellListMap.Box(unit_cell, dist_cutoff)
end
if isnothing(x0)
x = [ustrip.(box.unit_cell_max) .* rand(SVector{2, T}) for _ in 1:np]
else
x = x0
end
# Construct the cell list for the first time, to allocate 
cl = CellList(x, box; parallel=true, nbatches=number_of_batches)
return CellListMapNeighborFinder2D{2, T}(
nb_matrix, matrix_14, n_steps, dist_cutoff,
cl, CellListMap.AuxThreaded(cl), 
[NeighborList(0, [(0, 0, false)]) for _ in 1:CellListMap.nbatches(cl)])
end

function Molly.find_neighbors(s::System,
    nf::CellListMapNeighborFinder2D,
    current_neighbors=nothing,
    step_n::Integer=0;
    parallel::Bool=true)
!iszero(step_n % nf.n_steps) && return current_neighbors

if isnothing(current_neighbors)
neighbors = NeighborList()
else
neighbors = current_neighbors
end
aux = nf.aux
cl = nf.cl
neighbors.n = 0
neighbors_threaded = nf.neighbors_threaded
if parallel
for i in 1:length(neighbors_threaded)
neighbors_threaded[i].n = 0
end
else
neighbors_threaded[1].n = 0
end

box = CellListMap.Box(s.box_size, nf.dist_cutoff; lcell=1)
cl = UpdateCellList!(s.coords, box, cl, aux; parallel=parallel)

map_pairwise!(
(x, y, i, j, d2, pairs) -> push_pair!(pairs, i, j, nf.nb_matrix, nf.matrix_14),
neighbors, box, cl;
reduce=reduce_pairs,
output_threaded=neighbors_threaded,
parallel=parallel
)

nf.cl = cl
return neighbors
end