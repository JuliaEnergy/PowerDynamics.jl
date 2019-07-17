using LightGraphs: edges, nv, AbstractGraph
using NetworkDynamics: network_dynamics

"""
Encapsulates nodes & lines of the power grid and a graph connecting both.
"""
struct PowerGrid
    graph:: G where G <: AbstractGraph
    nodes
    lines
end

function PowerGrid(nodes, lines)
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, l.from, l.to) for l in lines]
    PowerGrid(graph, nodes, lines)
end

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), pg.nodes))
