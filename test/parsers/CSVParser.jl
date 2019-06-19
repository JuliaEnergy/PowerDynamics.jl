using Test: @test, @testset
using LightGraphs: edges, Edge


power_grid = read_network_from_csv(joinpath(@__DIR__, "test_nodes.csv"), joinpath(@__DIR__, "test_lines.csv"))

@test length(power_grid.nodes) == 3
@test length(power_grid.lines) == 1

@test collect(edges(power_grid.graph)) == [Edge(1,2)]
