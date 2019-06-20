using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: PVAlgebraic, construct_vertex, symbolsof

include("NodeTestBase.jl")

@testset "PVAlgebraic" begin
    @syms P real=true
    @syms V positive=true

    pv = PVAlgebraic(P=P, V=V)
    pv_vertex = construct_vertex(pv)
    @test symbolsof(pv) == [:u_r, :u_i]
    @test pv_vertex.mass_matrix == [0,0]

    smoketest_rhs(pv_vertex)
end
