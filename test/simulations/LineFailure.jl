using Test: @test
using LightGraphs: edges, Edge
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, LineFailure,LineFault, simulate, PowerGrid
using OrderedCollections: OrderedDict

Y = 0 + 5*im
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=1, to=2, Y=Y),StaticLine(from=2, to=3, Y=Y)]
grid = PowerGrid(nodes, lines)

linefault = LineFault(;from=1, to=2)
faulty_grid = linefault(grid)

@test length(faulty_grid.lines) == 1
@test collect(edges(faulty_grid.graph)) == [Edge(2,3)]

Y = 0 + 5*im
nodes = OrderedDict("bus1"=>SlackAlgebraic(U=1), "bus2"=>SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), "bus3"=>SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1))
lines = OrderedDict("line1"=>StaticLine(from="bus1", to="bus2", Y=Y),"line2"=>StaticLine(from="bus2", to="bus3", Y=Y))
grid = PowerGrid(nodes, lines)

linefailure = LineFailure(;line_name="line1")
faulty_grid = linefailure(grid)

@test length(faulty_grid.lines) == 1
@test collect(edges(faulty_grid.graph)) == [Edge(2,3)]