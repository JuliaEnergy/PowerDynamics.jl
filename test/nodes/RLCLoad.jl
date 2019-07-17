using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: RLCLoad, construct_vertex, symbolsof

include("NodeTestBase.jl")

@testset "RLCLoad" begin
    @syms R L C positive=true

    pv = RLCLoad(R=R, L=L, C=C)
    pv_vertex = construct_vertex(pv)
    @test symbolsof(pv) == [:u_r, :u_i, :u_C, :i_L]
    @test pv_vertex.mass_matrix == [0,0,1,1,1]

    smoketest_rhs(pv_vertex)
end
