#using Pkg
#Pkg.activate(".")
using PowerDynamics
using LightGraphs: SimpleGraph, add_edge!, add_vertices!, Edge
using MetaGraphs: MetaGraph, set_prop!


admittance = 1im/0.2 #+ 0.03


graph = MetaGraph(SimpleGraph())
node_list = []
append!(node_list, [SlackAlgebraic(U=1)])
append!(node_list, [RLCLoad(R=0.5,L=0.5,C=0.5)])
num_nodes = length(node_list)
add_vertices!(graph, num_nodes)
[set_prop!(graph, n, :node, node_list[n]) for n=1:num_nodes]

line_list =[]
line = StaticLine(Y=admittance)
push!(line_list, line)
add_edge!(graph, 1, 2);
set_prop!(graph, Edge(1, 2), :line, line)
powergrid = PowerGrid(graph, node_list, line_list)

operationpoint = find_operationpoint(powergrid)

result = simulate(Perturbation(1, :ω, Inc(0.2)), powergrid, operationpoint, timespan = (0.0,0.3))
plot_res(result, powergrid)
