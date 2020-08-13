using Test: @test, @testset
using PowerDynamics: SlackAlgebraic,StaticLine, SwingEqLVS, PowerGrid, systemsize, rhs
using LightGraphs: edges, Edge, nv, ne
using OrdinaryDiffEq: ODEFunction
using OrderedCollections: OrderedDict


nodes_dict = OrderedDict("bus1"=>SlackAlgebraic(U=1), "bus2"=> SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1))
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=1, to=2, Y=0*5im)]
lines_dict = OrderedDict("line1"=>StaticLine(from="bus1", to="bus2", Y=0*5im))

power_grid = PowerGrid(nodes, lines)
power_grid_fromdict = PowerGrid(nodes_dict,lines_dict)

@test systemsize(power_grid) == 5 # -> u_r_sl, u_i_sl, u_r_sw, u_i_sw, omega_sw
@test systemsize(power_grid_fromdict) == 5 # -> u_r_sl, u_i_sl, u_r_sw, u_i_sw, omega_sw

@test nv(power_grid.graph) == 2
@test nv(power_grid_fromdict.graph) == 2

@test ne(power_grid.graph) == 1
@test ne(power_grid_fromdict.graph) == 1

@test collect(edges(power_grid.graph)) == [Edge(1,2)]
@test collect(edges(power_grid_fromdict.graph)) == [Edge(1,2)]

@test rhs(power_grid) isa(ODEFunction)
@test rhs(power_grid_fromdict) isa(ODEFunction)
