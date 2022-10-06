using Test: @testset, @test
using PowerDynamics: FluctuationNode, construct_vertex, symbolsof
using LinearAlgebra: diag

include("../NodeTestBase.jl")

@testset "FluctuationNode Tests" begin
    p, q, ω = rand_real(3)

    fn = FluctuationNode(t -> p*sin(ω*t), t -> q)
    fn_vertex = construct_vertex(fn)

    @test symbolsof(fn) == [:u_r, :u_i]
    @test fn_vertex.mass_matrix |> diag == [0,0]

    smoketest_rhs(fn_vertex)
end