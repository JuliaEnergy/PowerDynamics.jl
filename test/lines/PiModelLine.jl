using Test: @testset, @test
using PowerDynamics: PiModelLine, construct_edge
using NetworkDynamics: StaticEdge

line = PiModelLine(from=1, to=2, y=10*im, y_shunt_km = 10, y_shunt_mk= 100)
edge = construct_edge(line)

@testset "PiModelLine construct_edge should return StaticEdge" begin
    @test isa(edge, StaticEdge)
end

@testset "PiModelLine should run edge function without errors" begin
        e = zeros(4)
        v_s = [10; 5*im]
        v_d = [12; 2*im]
        # assure function call does not explode!
        edge.f!(e, v_s, v_d, 0, 0)

        @test e[1] == 50
        @test e[2] == -50
        @test e[3] == 1000
        @test e[4] == 50
end
