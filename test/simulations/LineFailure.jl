using Test: @test
using LightGraphs: edges, Edge
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, Perturbation, LineFaailure, simulate

Y = 0 + 5*im
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(nodes, lines)

line_fault = LineFault(from=1, to=2)
faulty_grid = line_fault(grid)

@test length(faulty_grid.lines) == 1
@test collect(edges(faulty_grid.graph)) == [Edge(1,3)]
