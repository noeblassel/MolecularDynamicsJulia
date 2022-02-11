struct StateLogger
    n_steps::Int
    prefix::AbstractString
end

StateLogger(n_steps::Integer) = StateLogger(n_steps, "logfile")
StateLogger(n_steps::Integer, file_prefix::AbstractString) = StateLogger(n_steps, file_prefix)

function Molly.log_property!(logger::StateLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        save_reduced_lj_state(s,logger.prefix*"_$(step_n).txt")
    end
end


function save_reduced_lj_state(s::System,filename::AbstractString)
    file=open(filename,"w")
    println(file,"$(s.box_size[1]) $(s.box_size[2]) $(s.box_size[3])")
    for i = 1:length(s)
        x=s.coords[i]
        v=s.velocities[i]
        println(file,"$(x[1]) $(x[2]) $(x[3]) $(v[1]) $(v[2]) $(v[3])")
    end
    close(file)
end