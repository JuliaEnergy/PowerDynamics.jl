using PowerDynBase
using NetworkDynamics
using LightGraphs
using Test
include("../../src/solve/PowerGridSolutions.jl")
include("../../src/solve/solve.jl")


function tests()
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        nodes = [PQAlgebraic(S=1+im*1), SwingEq(H=1., P=-1, D=1, Ω=50), SwingEq(H=1., P=-1, D=1, Ω=50)]
        lines = [StaticLine(Y=1.0im) for e in edges(g)]
        nd = network_dynamics(map(construct_node_dynamics, nodes), map(construct_edge, lines), g)
        powergrid = PowerGrid(nd, g, nodes, lines)

        x0 = ones(systemsize(powergrid))
        timespan = (0.,10.)
        sol = solve(powergrid, x0, timespan)
end

tests()
