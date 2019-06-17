using DifferentialEquations: ODEProblem, ODEFunction
using MetaGraphs, LightGraphs

"""
Encapsulates nodes & lines of the power grid and a graph connecting both.
"""
struct PowerGrid
    graph:: G where G <: AbstractGraph
    nodes
    lines
end

function PowerGrid(graph::G) where G <: AbstractMetaGraph
    nodes = [get_prop(graph, n, :node) for n=1:nv(graph)]
    lines = [get_prop(graph, e, :line) for e in edges(graph)]
    PowerGrid(graph, nodes, lines)
end

function ode_function(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), pg.nodes))
