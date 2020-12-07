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
    line_array = collect(values(lines))
    line_keys = collect(keys(lines))
    
    # assert that keys are consistent
    @assert all([l.from ∈ bus_array && l.to ∈ bus_array for l in values(lines)])
    
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, findfirst(bus_array .== l.from), findfirst(bus_array .== l.to)) for (key,l) in lines]
    
    pg=PowerGrid(graph, nodes, lines)

    sources = [findfirst(bus_array .== l.from) for l in line_array]
    dest = [findfirst(bus_array .== l.to) for l in line_array]
    sorted_lines = OrderedDict()#deepcopy(lines)
    for (j,edge) in enumerate(collect(edges(pg.graph)))
        try 
            index = max(findfirst(sources.==edge.src),findfirst(dest.==edge.dst))
            sorted_lines[line_keys[index]]=lines[line_keys[index]]
        catch error_message
            index = max(findfirst(sources.==edge.dst),findfirst(dest.==edge.src))
            try  sorted_lines[line_keys[index]]=lines[line_keys[index]]
            catch
                println("no nodes matching the graph found")
            end
        end
    end

    PowerGrid(pg.graph,pg.nodes,sorted_lines)
end

function PowerGrid(nodes::Array, lines::Array)
    # assert that keys are consistent
    @assert all([l.from isa Int && 1 <= l.from <= length(nodes) for l in lines])
    @assert all([l.to isa Int && 1 <= l.from <= length(nodes) for l in lines])
  
    graph = SimpleGraph(length(nodes))
    [add_edge!(graph, l.from, l.to) for l in lines]
    pg=PowerGrid(graph, nodes, lines)

    bus_array=collect(keys(pg.nodes))
    lines = collect(values(pg.lines))
    sources = [findfirst(bus_array .== l.from) for l in lines]
    dest = [findfirst(bus_array .== l.to) for l in lines]
    sorted_lines = deepcopy(lines)
    for (j,edge) in enumerate(collect(edges(pg.graph)))
        try sorted_lines[j]=lines[max(findfirst(sources.==edge.src),findfirst(dest.==edge.dst))]
        catch error_message
            try  sorted_lines[j]=lines[max(findfirst(sources.==edge.dst),findfirst(dest.==edge.src))]
            catch
                println("no nodes matching the graph found")
            end
        end
    end

    PowerGrid(pg.graph,pg.nodes,sorted_lines)
end

function rhs(pg::PowerGrid)

    network_dynamics(map(construct_vertex, collect(values(pg.nodes))), map(construct_edge, collect(values(pg.lines))), pg.graph)
end

"""
Returns the total size of dynamical variables of the whole powergrid
"""
@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), collect(values(pg.nodes))))
