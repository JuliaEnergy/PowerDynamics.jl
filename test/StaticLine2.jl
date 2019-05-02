using Test
using PowerDynBase
using NetworkDynamics

function tests()

    edge::StaticEdge = StaticLine2!(Y=10*im)

    @testset "StaticLine2! should convert to StaticEdge" begin
        @test isa(edge, StaticEdge)
    end

    @testset "StaticLine2! should run edge function without errors" begin
        e = []
        v_s = [10; 5*im]
        v_d = [12; 2*im]
        # assure function call does not explode!
        edge.f!(e, v_s, v_d, 0, 0)
    end
end

tests()
