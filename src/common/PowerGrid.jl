using LightGraphs: src, dst, edges, nv, AbstractGraph
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
    bus_keys=collect(keys(nodes))
    line_array = collect(values(lines))
    line_keys = collect(keys(lines))
    
    # assert that keys are consistent
    @assert all([l.from ∈ bus_keys && l.to ∈ bus_keys for l in values(lines)]) "invalid node key given in lines"

    # Assert ordering of from/to according to index in bus_keys to assure
    # to assure compatibility with LighGraphs.
    @assert all([findfirst(bus_keys .== l.from) < findfirst(bus_keys .== l.to) for l in values(lines)]) "the pairs (from, to) need to be ordered according to the index of the corresponding keys in the node dict"

    graph = SimpleGraph(length(nodes))
    for (key,l) in lines
        add_edge!(graph, findfirst(bus_keys .== l.from), findfirst(bus_keys .== l.to)) 
    end

    sorted_lines = OrderedDict()#deepcopy(lines)
    for e in edges(graph)
        original_key = findfirst( x ->  bus_keys[src(e)] == x.from && bus_keys[dst(e)] == x.to, lines)
        # We can maintain the same keys but simply reorder the OrderedDict since 
        # the entries appear in the order they're added.
        sorted_lines[original_key] = lines[original_key]
    end

    PowerGrid(graph,nodes,sorted_lines)
end

function PowerGrid(nodes::Array, lines::Array)
    # assert that keys are consistent
    @assert all([l.from isa Int for l in lines]) "`from` should be of type `Int`"
    @assert all([l.to isa Int for l in lines]) "`to` should be of type `Int`"
    @assert all([1 <= l.from <= length(nodes) for l in lines]) "numerical index needs to be between 1 and the number of nodes"
    @assert all([1 <= l.to <= length(nodes) for l in lines]) "numerical index needs to be between 1 and the number of nodes"

    # We should enforce ordering of from/to to comply with Lightgraphs.jl. 
    # This could otherwise lead to problems for unsymmetric line types.
    @assert all([l.from < l.to for l in lines]) "the pairs (from, to) need to be ordered according to index value"
  
    graph = SimpleGraph(length(nodes))
    for l in lines
        add_edge!(graph, l.from, l.to)
    end
    
    # We sort the lines to be in accordance with the indexing in the graph.
    sorted_lines = similar(lines)
    for (idx, e) in enumerate(edges(graph))
        original_idx = findfirst( x ->  src(e) == x.from && dst(e) == x.to, lines)
        sorted_lines[idx] = lines[original_idx]
    end

    PowerGrid(graph, nodes, sorted_lines)

    # pg=PowerGrid(graph, nodes, lines)

    # bus_array=collect(keys(pg.nodes))
    # lines = collect(values(pg.lines))
    # sources = [findfirst(bus_array .== l.from) for l in lines]
    # dest = [findfirst(bus_array .== l.to) for l in lines]
    # sorted_lines = deepcopy(lines)
    # for (j,edge) in enumerate(collect(edges(pg.graph)))
    #     try sorted_lines[j]=lines[max(findfirst(sources.==edge.src),findfirst(dest.==edge.dst))]
    #     catch error_message
    #         try  sorted_lines[j]=lines[max(findfirst(sources.==edge.dst),findfirst(dest.==edge.src))]
    #         catch
    #             println("no nodes matching the graph found")
    #         end
    #     end
    # end

    # PowerGrid(pg.graph,pg.nodes,sorted_lines)
end

function rhs(pg::PowerGrid)

    network_dynamics(map(construct_vertex, collect(values(pg.nodes))), map(construct_edge, collect(values(pg.lines))), pg.graph)
end

"""
Returns the total size of dynamical variables of the whole powergrid
"""
@views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), collect(values(pg.nodes))))
