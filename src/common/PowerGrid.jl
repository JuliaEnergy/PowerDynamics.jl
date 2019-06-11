using DifferentialEquations: ODEProblem, ODEFunction
using MetaGraphs, LightGraphs

struct PowerGrid
    graph:: G where G <: AbstractMetaGraph
    nodes
    lines
end

function PowerGrid(graph::G) where G <: AbstractMetaGraph
    nodes = [get_prop(graph, n, :node) for n=1:nv(graph)]
    lines = [get_prop(graph, e, :line) for e in edges(graph)]
    PowerGrid(graph, nodes, lines)
end

function ode_function(pg::PowerGrid)
    network_dynamics(map(construct_node_dynamics, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), pg.nodes))
