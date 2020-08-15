using Test: @test, @testset, @test_throws
using PowerDynamics: SlackAlgebraic,StaticLine, SwingEqLVS, PowerGrid, systemsize, rhs
using LightGraphs: edges, Edge, nv, ne
using OrdinaryDiffEq: ODEFunction
using OrderedCollections: OrderedDict

##

nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
dnodes = ["bus$i" => nodes[i] for i in 1:length(nodes)]
lines = [StaticLine(from=1, to=2, Y=0*5im)]
dlines = ["line$i" => lines[i] for i in 1:length(lines)]
lines_string = [StaticLine(from="bus1", to="bus2", Y=0*5im)]
dlines_string = ["line$i" => lines_string[i] for i in 1:length(lines_string)]

##

# test array constructor
power_grid = PowerGrid(nodes, lines)
@test power_grid isa PowerGrid

@test_throws AssertionError PowerGrid(nodes, lines_string)

##

# test dict constructor
power_grid_fromdict = PowerGrid(OrderedDict(dnodes), OrderedDict(dlines_string))
@test power_grid_fromdict isa PowerGrid

@test_throws AssertionError PowerGrid(OrderedDict(dnodes), OrderedDict(dlines))

## test invalid input

@test_throws ErrorException PowerGrid(Dict(dnodes), Dict(dlines_string))
@test_throws ErrorException PowerGrid(OrderedDict(dnodes), lines_string)
@test_throws ErrorException PowerGrid(nodes, OrderedDict(dlines))

##
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
