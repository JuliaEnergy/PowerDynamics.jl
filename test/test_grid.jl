using PowerDynBase: SlackAlgebraic, SwingEqLVS, PowerGrid, StaticLine
using LightGraphs: SimpleGraph, add_edge!, edges


function create_grid()
    graph = SimpleGraph(3)
    nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
    add_edge!(graph, 1, 2)
    add_edge!(graph, 1, 3)
    lines = [StaticLine(Y=0.0 - 5.0im) for e in edges(graph)]
    PowerGrid(graph, nodes, lines)
end
