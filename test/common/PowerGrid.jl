using Test: @test, @testset
using PowerDynamics: SlackAlgebraic, SwingEqLVS, PowerGrid
using LightGraphs: edges, Edge, nv, ne
using OrdinaryDiffEq: ODEFunction


nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=1, to=2, Y=0*5im)]

power_grid = PowerGrid(nodes, lines)

@test systemsize(power_grid) == 5 # -> u_r_sl, u_i_sl, u_r_sw, u_i_sw, omega_sw
@test nv(power_grid.graph) == 2
@test ne(power_grid.graph) == 1
@test collect(edges(power_grid.graph)) == [Edge(1,2)]

@test rhs(power_grid) isa(ODEFunction)
