using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: WindTurbineGenType4_RotorControl, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "WindGenType 4 with Rotor Control Tests" begin
    #T_L,T_H,K_P,ΔP_max,K_PLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_g1,K_g2,K_r1,K_r2
    #[θ_PLL,dθ_PLL],[e_Idθ,de_Idθ],[ω,dω],[e_IP,de_IP],[z,dz],[e_IV,de_IV],[u_dc,du_dc],[i_q,di_q],[u_tref,du_tref],[ω_r,dω_r]
    @syms T_L T_H K_P ΔP_max K_PLL C J  u_dcref K_Q K_v K_g1 K_g2 K_r1 K_r2 positive=true
    @syms Q_ref ω_rref P real=true
    @syms θ_PLL dθ_PLL e_Idθ de_Idθ ω dω  e_IP de_IP z dz e_IV de_IV u_dc du_dc i_q di_q u_tref du_tref ω_r dω_r real=true
    @syms k_P_invalid k_PLL_INVALID negative=true

    @test_throws AssertionError construct_vertex(WindTurbineGenType4_RotorControl(T_L=T_L,T_H=T_H,K_P=K_P,ΔP_max=ΔP_max,K_PLL=k_PLL_INVALID,C=C,J=J,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2))
    @test_throws AssertionError construct_vertex(WindTurbineGenType4_RotorControl(T_L=T_L,T_H=T_H,K_P=k_P_invalid,ΔP_max=ΔP_max,K_PLL=K_PLL,C=C,J=J,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2))

    wind_gen = WindTurbineGenType4_RotorControl(T_L=T_L,T_H=T_H,K_P=K_P,ΔP_max=ΔP_max,K_PLL=K_PLL,C=C,J=J,u_dcref=u_dcref,K_Q=K_Q,K_v=K_v,K_g1=K_g1,K_g2=K_g2,K_r1=K_r1,K_r2=K_r2)
    wind_gen_vertex = construct_vertex(wind_gen)

    @test symbolsof(wind_gen) == [:u_r,:u_i,:θ_PLL,:e_Idθ,:ω,:e_IP,:z,:e_IV,:u_dc,:i_q,:u_tref,:ω_r]
    @test wind_gen_vertex.mass_matrix == [0,0,1,1,1,1,1,1,1,1,1,1]

    smoketest_rhs(wind_gen_vertex, int_x=[θ_PLL,e_Idθ,e_Idθ,ω,e_IP,z,e_IV,u_dc,i_q,u_tref,ω_r], int_dx=[dθ_PLL,de_Idθ,de_Idθ,dω,de_IP,dz,de_IV,du_dc,di_q,du_tref,dω_r])
end
