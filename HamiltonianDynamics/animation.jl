using GLMakie,Makie
using Molly

function visualize_alt(coord_logger,
    box_size,
    out_filepath::AbstractString;
    framerate::Integer=30,
    color=:purple,
    markersize=20.0,
    linewidth=2.0,
    transparency=true,
    kwargs...)
coords_start = first(coord_logger.coords)
dims = length(first(coords_start))
if dims == 3
PointType = Point3f
elseif dims == 2
PointType = Point2f
else
throw(ArgumentError("Found $dims dimensions but can only visualize 2 or 3 dimensions"))
end

scene = Scene()
positions = Observable(PointType.(ustrip_vec.(coords_start)))
scatter!(scene, positions; color=color, markersize=markersize,
transparency=transparency, kwargs...)

dist_unit = unit(first(first(coords_start)))
box_size_conv = ustrip.(dist_unit, box_size)
println(box_size_conv)
xlims!(scene, 0.0, box_size_conv[1])
ylims!(scene, 0.0, box_size_conv[2])
dims == 3 && zlims!(scene, 0.0, box_size_conv[3])

GLMakie.record(scene, out_filepath, eachindex(coord_logger.coords); framerate=framerate) do frame_i
coords = coord_logger.coords[frame_i]
end
end
