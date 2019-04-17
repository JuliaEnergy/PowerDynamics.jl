using Test
using PowerDynBase
using NetworkDynamics

@testset "StaticLine2! should convert to StaticEdge" begin
    edge::StaticEdge = StaticLine2!(Y=10*im)
    @test isa(edge, StaticEdge)
end
