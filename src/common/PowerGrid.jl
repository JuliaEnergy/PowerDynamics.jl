using LightGraphs: edges, nv, AbstractGraph
using NetworkDynamics: network_dynamics

"""
```Julia
Powergrid(graph, nodes, lines)
```

Our model describes a powergrid as an undirected graph consisting
of edges that represent electrical lines and vertices that represent
specific electrical nodes i.e. generators, inverters etc.

"""
struct PowerGrid
    graph:: G where G <: AbstractGraph
    nodes
    lines
end

"""
```Julia
Powergrid(nodes, lines)
```

creates a [`PowerGrid`](@ref) from nodes and lines. The underlying graph
is created automatically.

"""
function PowerGrid(nodes::Dict, lines::Dict)
    graph = SimpleGraph(length(nodes))
    bus_array=collect(keys(nodes))
    [add_edge!(graph, findfirst(x->x==l.from, bus_array), findfirst(x->x==l.to, bus_array)) for (key,l) in lines]
    PowerGrid(graph, nodes, lines)
end

function PowerGrid(nodes::Array, lines::Array)
  
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, l.from, l.to) for l in lines]
    PowerGrid(graph, nodes, lines)
end

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, collect(values(pg.nodes))), map(construct_edge, collect(values(pg.lines))), pg.graph)
end

"""
Returns the total size of dynamical variables of the whole powergrid
"""
@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), collect(values(pg.nodes))))
