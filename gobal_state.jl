module GlobalState

export velocity_queue, initialize_velocity_queue, update_velocity_queue

global velocity_queue = []

function initialize_velocity_queue(size::Int)
    global velocity_queue = fill(2, size)
end

function update_velocity_queue(values::Vector{Float64})
    global velocity_queue = values
end

end # module
