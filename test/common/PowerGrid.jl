using Test: @test, @testset
using PowerDynamics: SlackAlgebraic, SwingEqLVS, PowerGrid
using LightGraphs: SimpleGraph, add_edge!, edges, Edge
using MetaGraphs: MetaGraph, set_prop!
using DifferentialEquations: ODEFunction


nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
graph = MetaGraph(SimpleGraph(2))
[set_prop!(graph, n, :node, nodes[n]) for n=1:length(nodes)]

add_edge!(graph, 1, 2);
line = StaticLine(Y=0*5im)
set_prop!(graph, Edge(1, 2), :line, line)

power_grid = PowerGrid(graph)

@test systemsize(power_grid) == 5 # -> u_r_sl, u_i_sl, u_r_sw, u_i_sw, omega_sw
@test power_grid.nodes == nodes
@test power_grid.lines == [line]

@test ode_function(power_grid) isa(ODEFunction)
