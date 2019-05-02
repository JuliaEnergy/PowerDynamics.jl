using Test
using PowerDynBase
using NetworkDynamics

function tests()

    line = StaticLine(Y=10*im)
    edge = construct_edge(line)

    @testset "construct_edge should return StaticEdge" begin
        @test isa(edge, StaticEdge)

        e = []
        v_s = [10; 5*im]
        v_d = [12; 2*im]
        edge.f!(e, v_s, v_d, 0, 0)
    end

    @testset "should run edge function without errors" begin
        e = []
        v_s = [10; 5*im]
        v_d = [12; 2*im]
        # assure function call does not explode!
        edge.f!(e, v_s, v_d, 0, 0)
    end
end

tests()
