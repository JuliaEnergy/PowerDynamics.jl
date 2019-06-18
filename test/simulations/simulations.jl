using Test: @test
using LightGraphs: SimpleGraph, add_edge!, edges, Edge
using MetaGraphs: MetaGraph, set_prop!

Y = 0 + 5*im
graph = SimpleGraph(3)
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
add_edge!(graph, 1, 2)
lines = [StaticLine(Y=Y) for e in edges(graph)]
grid = PowerGrid(graph, nodes, lines)
@syms u_Sl u_Sw1
omega1 = 0.01
omega1_inc = 0.2
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1])
state2 = Perturbation(2, :ω, Inc(omega1_inc))(state)

@test state2[2, :ω] == omega1 + omega1_inc
