using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: RLCLoad, construct_vertex, symbolsof

include("NodeTestBase.jl")

@testset "RLCLoad" begin
    @syms R L C positive=true
    @syms u_C i_L omega domega du_C di_L real=true

    pv = RLCLoad(R=R, L=L, C=C)
    pv_vertex = construct_vertex(pv)
    @test symbolsof(pv) == [:u_r, :u_i, :u_C, :i_L, :ω]
    @test pv_vertex.mass_matrix == [0,0,1,1,1]

    smoketest_rhs(pv_vertex, int_x=[u_C, i_L, omega], int_dx=[domega, du_C, di_L])
end
