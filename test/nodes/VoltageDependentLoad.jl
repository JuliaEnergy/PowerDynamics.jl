using Test: @testset, @test
using PowerDynamics: PQAlgebraic, SlackAlgebraic, StaticLine, VoltageDependentLoad, construct_vertex, symbolsof, find_operationpoint, PowerGrid

include("NodeTestBase.jl")

@testset "VoltageDependentLoad" begin
    P, Q = rand_real(2)
    U = 1
    A = 0.3
    B = 0.3
    load = VoltageDependentLoad(P = P, Q = Q, U = U, A = A, B = B)
    load_vertex = construct_vertex(load)
    @test symbolsof(load) == [:u_r, :u_i]
    @test load_vertex.mass_matrix == [0,0]

    smoketest_rhs(load_vertex)
end

@testset "PQAlgebraic as a corner case of VoltageDependentLoad" begin
    P, Q = rand_real(2)
    U = 1
    A = 0.
    B = 0.

    slack = SlackAlgebraic(U = U)
    line = StaticLine(;from = 1, to = 2, Y = 10 / (0.1152 + im * 0.0458))

    load = VoltageDependentLoad(P = P, Q = Q, U = U, A = A, B = B)
    pq = PQAlgebraic(P = P, Q = Q)

    pg_pq = PowerGrid([slack, pq], [line,])
    pg_load = PowerGrid([slack, load], [line,])

    op_load = find_operationpoint(pg_load)
    op_pq = find_operationpoint(pg_pq)

    # in the steady state, we expect equal currents/voltage
    @test op_pq[2, :i] == op_load[2, :i]
    @test op_pq[2, :v] == op_load[2, :v]

end
