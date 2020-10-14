using Test: @testset, @test
using PowerDynamics: SlackAlgebraic, construct_vertex, symbolsof
using LinearAlgebra: diag

include("NodeTestBase.jl")

@testset "SlackAlgebraic" begin
    U = rand(ComplexF64)
    slack = SlackAlgebraic(U=U)
    slack_vertex = construct_vertex(slack)

    @test symbolsof(slack) == [:u_r, :u_i]
    @test slack_vertex.mass_matrix |> diag == [0,0]

    smoketest_rhs(slack_vertex)
end
