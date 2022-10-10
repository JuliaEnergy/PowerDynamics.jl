using Test: @testset, @test
using PowerDynamics: PiModelLine, Transformer, PiModel, construct_edge
using NetworkDynamics: StaticEdge

trafo = Transformer(from=1, to=2, Y=10*im, T_ratio=2.)
edge = construct_edge(trafo)

@testset "Transformer construct_edge should return StaticEdge" begin
    @test isa(edge, StaticEdge)
    @test edge.dim == 4
end

@testset "Transformer should run edge function without errors" begin
        e = zeros(4)
        v_s = [1.; 0]
        v_d = [sqrt(2)/2; sqrt(2)/2]
        # assure function call does not explode!
        edge.f(e, v_s, v_d, 0, 0)

        I_km = -abs2(trafo.T_ratio) * trafo.Y * complex(v_s...) + conj(trafo.T_ratio) * trafo.Y * complex(v_d...)
        I_mk = trafo.Y * complex(v_d...) - trafo.T_ratio * trafo.Y * complex(v_s...)

        @test e[3] ≈ -real(I_km)
        @test e[4] ≈ -imag(I_km)
        @test e[1] ≈ real(I_mk)
        @test e[2] ≈ imag(I_mk)
end

@testset "Transformer should return the same result as PiModelLine" begin
        trafo = Transformer(;from=1, to=2, Y=1/(0.1152 + im*0.0458), T_ratio=1.)
        pml = PiModelLine(;from=1, to=2, Y = 1/(0.1152 + im*0.0458), Y_shunt_km = 0.,  Y_shunt_mk = 0.)

        arr = zeros(4)
        ed = construct_edge(trafo)
        ed.f(arr, [0.1,0.2], [0.3,0.4], 0, 0)

        arr2 = zeros(4)
        ed2 = construct_edge(pml)
        ed2.f(arr2, [0.1,0.2], [0.3,0.4], 0, 0)

        @test isapprox(arr,arr2)
end
