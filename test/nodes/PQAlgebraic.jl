using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: PQAlgebraic, construct_vertex, symbolsof

include("NodeTestBase.jl")

@testset "PQAlgebraic" begin
    @syms S
    pq = PQAlgebraic(S=S)
    pq_vertex = construct_vertex(pq)
    @test symbolsof(pq) == [:u_r, :u_i]
    @test pq_vertex.mass_matrix == [0,0]

    smoketest_rhs(pq_vertex)
end
