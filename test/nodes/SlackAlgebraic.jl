using Test
using SymPy
using PowerDynBase

include("NodeTestBase.jl")

@testset "SlackAlgebraic" begin
    @syms U
    slack = SlackAlgebraic(U=U)
    slack_vertex = construct_vertex(slack)

    @test symbolsof(slack) == [:u_r, :u_i]
    @test slack_vertex.mass_matrix == [0,0]

    smoketest_rhs(slack_vertex)
end
