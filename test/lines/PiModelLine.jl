using Test: @testset, @test
using PowerDynamics: PiModelLine, StaticLine, PiModel, construct_edge
using NetworkDynamics: StaticEdge

line = PiModelLine(from=1, to=2, y=10*im, y_shunt_km = 10, y_shunt_mk= 100)
edge = construct_edge(line)

@testset "PiModelLine construct_edge should return StaticEdge" begin
    @test isa(edge, StaticEdge)
end

@testset "PiModelLine should run edge function without errors" begin
        e = zeros(4)
        v_s = [10; 5]
        v_d = [12; 2]
        # assure function call does not explode!
        edge.f(e, v_s, v_d, 0, 0)

        @test e[1] == 1230
        @test e[2] == 220
        @test e[3] == 70
        @test e[4] == 30
end

@testset "PiModelLine should return the same result as StaticLine" begin
        sl = StaticLine(;from=1, to=4, Y = 1/(0.1152 + im*0.0458))
        pml = PiModelLine(;from=1, to=4, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.)

        arr = zeros(4)
        ed = construct_edge(sl)
        ed.f(arr, [0.1,0.2], [0.3,0.4], 0, 0)

        arr2 = zeros(4)
        ed2 = construct_edge(pml)
        ed2.f(arr2, [0.1,0.2], [0.3,0.4], 0, 0)

        @test isapprox(arr,arr2)
end
