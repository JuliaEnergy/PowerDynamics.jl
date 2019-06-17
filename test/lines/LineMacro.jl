using Test
using PowerDynBase
using NetworkDynamics

line = StaticLine(Y=10*im)
edge = construct_edge(line)

@testset "LineMacro construct_edge should return StaticEdge" begin
    @test isa(edge, StaticEdge)
end

@testset "LineMacro should run edge function without errors" begin
        e = zeros(4)
        v_s = [10; 5*im]
        v_d = [12; 2*im]
        # assure function call does not explode!
        edge.f!(e, v_s, v_d, 0, 0)

        @test e[1] == 0
        @test e[2] == 50
        @test e[3] == 0
        @test e[4] == 50
end
