using Test: @testset, @test
using PowerDynamics: RLCLoad, construct_vertex, symbolsof
using LinearAlgebra: diag

include("../NodeTestBase.jl")


@testset "RLCLoad" begin
    R,L,C = rand_positive(3)
    u_Cr, u_Ci,i_Lr,i_Li,du_Cr,du_Ci,di_Lr,di_Li,omega,domega  = rand_real(10)

    pv = RLCLoad(R=R, L=L, C=C)
    pv_vertex = construct_vertex(pv)
    @test symbolsof(pv) == [:u_r, :u_i, :u_Cr, :u_Ci, :i_Lr,:i_Li, :ω]
    @test pv_vertex.mass_matrix |> diag == [0,0,1,1,1,1,1]

    smoketest_rhs(pv_vertex, int_x=[u_Cr, u_Ci, i_Lr,i_Li, omega], int_dx=[domega, du_Cr,du_Ci, di_Lr,di_Li])
end
