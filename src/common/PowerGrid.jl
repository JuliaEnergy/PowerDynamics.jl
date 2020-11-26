using LightGraphs: edges, nv, AbstractGraph
using NetworkDynamics: network_dynamics
using OrderedCollections: OrderedDict

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

creates a [`PowerGrid`](@ref) from nodes and lines (either given as a list or as a dictionay). 
The underlying graph is created automatically.

"""
function PowerGrid(nodes, lines)
    throw(error("Please supply both bus and line components either in an `Array` or `OrderedDict` (e.g. from the package `OrderedCollections`)."))
end

function PowerGrid(nodes::OrderedDict, lines::OrderedDict)
    bus_array=collect(keys(nodes))
    
    # assert that keys are consistent
    @assert all([l.from ∈ bus_array && l.to ∈ bus_array for l in values(lines)])
    
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, findfirst(bus_array .== l.from), findfirst(bus_array .== l.to)) for (key,l) in lines]
    PowerGrid(graph, nodes, lines)
end

function PowerGrid(nodes::Array, lines::Array)
    # assert that keys are consistent
    @assert all([l.from isa Int && 1 <= l.from <= length(nodes) for l in lines])
    @assert all([l.to isa Int && 1 <= l.from <= length(nodes) for l in lines])
  
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, l.from, l.to) for l in lines]
    PowerGrid(graph, nodes, lines)
end

function rhs(pg::PowerGrid)
    sorted_lines = collect(values(pg.lines))
    network_dynamics(map(construct_vertex, collect(values(pg.nodes))), map(construct_edge, sorted_lines, pg.graph)
end

"""
Returns the total size of dynamical variables of the whole powergrid
"""
@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), collect(values(pg.nodes))))
