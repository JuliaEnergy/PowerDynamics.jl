using Test: @testset, @test
using PowerDynamics: PQAlgebraic, construct_vertex, symbolsof
using LinearAlgebra: diag

include("NodeTestBase.jl")

@testset "PQAlgebraic" begin
    P,Q = rand_real(2)
    pq = PQAlgebraic(P=P,Q=Q)
    pq_vertex = construct_vertex(pq)
    @test symbolsof(pq) == [:u_r, :u_i]
    @test pq_vertex.mass_matrix |> diag == [0,0]

    smoketest_rhs(pq_vertex)
end
