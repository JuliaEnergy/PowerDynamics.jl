using DifferentialEquations: ODEProblem, ODEFunction
using LightGraphs

struct PowerGrid
    network_dynamics
    graph:: AbstractGraph
    nodes:: Array{Any} # TODO abstract base type here?
    lines:: Array{Any} # TODO abstract base type here?
end

@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), pg.nodes))
