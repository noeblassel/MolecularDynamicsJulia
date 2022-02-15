function calculate_simulation(s::NBodySimulation, alg_type::SymplecticEuler, args...; kwargs...)
    cb = obtain_callbacks_for_so_ode_problem(s)
    solution = solve(SecondOrderODEProblem(s), alg_type, args...; callback=cb, save_everystep = isempty(cb), kwargs...)
    return SimulationResult(solution, s)
end