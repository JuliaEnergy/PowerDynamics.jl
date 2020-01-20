# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

"Abstract super type of all PowerDynamics.jl Errors."
abstract type PowerDynamicsError <: Exception end

"Error to be thrown if something goes wrong during the node dynamics construction."
struct NodeDynamicsError <: PowerDynamicsError
    msg::String
end

"Error to be thrown if something goes wrong when creating or modifying states."
struct StateError <: PowerDynamicsError
    msg::String
end

"Error to be thrown if something goes wrong during the operation point search."
struct OperationPointError <: PowerDynamicsError
    msg::String
end

"Error to be thrown if something goes wrong during power grid solving"
struct GridSolutionError <: PowerDynamicsError
    msg::String
end
