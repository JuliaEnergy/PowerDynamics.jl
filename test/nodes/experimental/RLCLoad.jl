using Test: @testset, @test
using PowerDynamics: RLCLoad, construct_vertex, symbolsof

include("../NodeTestBase.jl")


@testset "RLCLoad" begin
    R,L,C = rand_positive(3)
    u_C,i_L,omega,domega,du_C,di_L = rand_real(6)

    pv = RLCLoad(R=R, L=L, C=C)
    pv_vertex = construct_vertex(pv)
    @test symbolsof(pv) == [:u_r, :u_i, :u_C, :i_L, :ω]
    @test pv_vertex.mass_matrix == [0,0,1,1,1]

    smoketest_rhs(pv_vertex, int_x=[u_C, i_L, omega], int_dx=[domega, du_C, di_L])
end
