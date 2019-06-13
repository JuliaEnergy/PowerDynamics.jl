using Test: @testset, @test
using SymPy: @syms
using PowerDynBase: ExponentialRecoveryLoad, construct_vertex, symbolsof

@testset "ExponentialRecoveryLoad" begin
    @syms V0 Nps Npt Nqs Nqt Tp Tq positive=true
    @syms P0 Q0 Pd Qd real=true
    @syms x_p dx_p x_q dx_q real=true

    exponentialRecoveryLoad = ExponentialRecoveryLoad(P0=P0, Q0=Q0, Nps=Nps, Npt=Npt, Nqs=Nqs, Nqt=Nqt, Tp=Tp, Tq=Tq, V0=V0)
    exprec_vertex= construct_vertex(exponentialRecoveryLoad)

    @test symbolsof(exponentialRecoveryLoad) == [:u_r, :u_i, :x_p, :x_q]
    @test exprec_vertex.mass_matrix == [0,0,1,1]

    smoketest_rhs(exprec_vertex, int_x=[x_p,x_q], int_dx=[dx_p, dx_q])
end
