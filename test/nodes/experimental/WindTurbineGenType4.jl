using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: WindTurbineGenType4, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "WindGenType 4 Tests" begin
    @syms K_PLL C J P ω_rref K_Q K_v K_g1 K_g2 K_r1 K_r2 positive=true
    @syms Q_ref P u_dcref real=true
    @syms θ dθ e_IP de_IP e_IV de_IV u_dc du_dc i_q di_q u_tref du_tref ω_r dω_r real=true
    @syms C_INVALID k_PLL_INVALID negative=true

    @test_throws AssertionError construct_vertex(WindTurbineGenType4(K_PLL=k_PLL_INVALID,Q_ref=Q_ref,C=C,J=J,P=P,ω_rref=ω_rref,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2))
    @test_throws AssertionError construct_vertex(WindTurbineGenType4(K_PLL=K_PLL,Q_ref=Q_ref,C=C_INVALID,J=J,P=P,ω_rref=ω_rref,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2))

    wind_gen = WindTurbineGenType4(K_PLL=K_PLL,Q_ref=Q_ref,C=C,J=J,P=P,ω_rref=ω_rref,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2)
    wind_gen_vertex = construct_vertex(wind_gen)

    @test symbolsof(wind_gen) == [:u_r,:u_i,:θ,:e_IP,:e_IV,:u_dc,:i_q,:u_tref,:ω_r]
    @test wind_gen_vertex.mass_matrix == [0,0,1,1,1,1,1,1,1]

    smoketest_rhs(wind_gen_vertex, int_x=[θ,e_IP,e_IV,u_dc,i_q,u_tref, ω_r], int_dx=[dθ,de_IP,de_IV,du_dc,di_q,du_tref, dω_r])
end
