using Test: @test
using LightGraphs: SimpleGraph, add_edge!, edges, Edge
using MetaGraphs: MetaGraph, set_prop!

Y = 0 + 5*im
graph = MetaGraph(SimpleGraph(3))
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
[set_prop!(graph, n, :node, nodes[n]) for n=1:length(nodes)]

add_edge!(graph, 1, 2)
add_edge!(graph, 1, 3)
lines = [StaticLine(Y=Y) for e in edges(graph)]
set_prop!(graph, Edge(1, 2), :line, lines[1])
set_prop!(graph, Edge(1, 3), :line, lines[2])
grid = PowerGrid(graph, nodes, lines)
@syms u_Sl u_Sw1 u_Sw2
omega1 = 0.01
omega1_delta = 0.2
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])

state2 = Perturbation(2, :ω, Inc(omega1_delta))(state)
@test state2[2, :ω] == omega1 + omega1_delta

state3 = Perturbation(2, :ω, Dec(omega1_delta))(state)
@test state3[2, :ω] == omega1 - omega1_delta

line_fault = LineFault(from=1, to=2)
faulty_grid = line_fault(grid)

@test length(faulty_grid.lines) == 1
@test collect(edges(faulty_grid.graph)) == [Edge(1,3)]
