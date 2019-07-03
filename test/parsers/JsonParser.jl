using PowerDynamics: read_network_from_json
using Test: @test

grid = read_network_from_json(joinpath(@__DIR__, "grid.json"))
@test grid != nothing
@test grid.lines[1].from == 1
@test grid.lines[1].to == 2
@test grid.lines[1].Y == 0.1+5im
