using Test: @testset, @test
using PowerDynamics: ExponentialRecoveryLoad, construct_vertex, symbolsof
using LinearAlgebra: diag

include("NodeTestBase.jl")

@testset "ExponentialRecoveryLoad" begin
    V0,Nps,Npt,Nqs,Nqt,Tp,Tq = rand_positive(7)
    P0,Q0,Pd,Qd = rand_real(4)
    x_p,dx_p,x_q,dx_q = rand_real(4)

    exponentialRecoveryLoad = ExponentialRecoveryLoad(P0=P0, Q0=Q0, Nps=Nps, Npt=Npt, Nqs=Nqs, Nqt=Nqt, Tp=Tp, Tq=Tq, V0=V0)
    exprec_vertex= construct_vertex(exponentialRecoveryLoad)

    @test symbolsof(exponentialRecoveryLoad) == [:u_r, :u_i, :x_p, :x_q]
    @test exprec_vertex.mass_matrix |> diag == [0,0,1,1]

    smoketest_rhs(exprec_vertex, int_x=[x_p,x_q], int_dx=[dx_p, dx_q])
end
