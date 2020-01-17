using Test: @testset, @test
using PowerDynamics: StepTransformer, construct_edge
using NetworkDynamics: StaticEdge

T = StepTransformer(Uos=110.,
                Uus=20.,
                k=0,
                ssp=10.,
                stufe_us=0.625,
                Sr=25.,
                uk=12.,
                Pvk=25.,
                Pvl=0.)

edge = construct_edge(T)

@testset "LineMacro construct_edge should return StaticEdge" begin
    @test isa(edge, StaticEdge)
end

@testset "LineMacro should run edge function without errors" begin
        e = zeros(4)
        v_s = [110E3; 0.]
        v_d = [20.279E3; 2.48078E3 * im]
        # assure function call does not explode!
        edge.f!(e, v_s, v_d, 0, 0)

        @test e[1] == 0
        @test e[2] == 50
        @test e[3] == 0
        @test e[4] == 50
end
