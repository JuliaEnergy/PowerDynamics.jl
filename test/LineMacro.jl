using Test
using PowerDynBase
using NetworkDynamics

@testset "construct_edge should return StaticEdge" begin
    line = StaticLine(Y=10*im)
    edge = construct_edge(line)
    @test isa(edge, StaticEdge)
end
