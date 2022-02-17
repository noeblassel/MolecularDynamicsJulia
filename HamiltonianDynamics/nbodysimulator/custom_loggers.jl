using NBodySimulator

function virial(res::SimulationResult,t::Real)
    N=res.
end

function pressure(res::SimulationResult,t::Real)

    V=nothing

    if isa(res.simulation.boundary_conditions,CubicPeriodicBoundaryConditions)
        V=res.simulation.boundary_conditions.L^3
    elseif isa(res.simulation.boundary_conditions,PeriodicBoundaryConditions)
        Lxm,LxM,Lym,LyM,Lzm,LzM=res.simulation.boundary_conditions.boundary
        V=(LxM-Lxm)*(LyM-Lym)*(LzM-Lzm)
    else
        error("Pressure only supported in finite domains.")
    end

    
end