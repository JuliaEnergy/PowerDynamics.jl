using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: PVInverterWithFrequencyControl, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "PVInverter Tests" begin
    @syms k_PLL f f_s T_m k_P positive=true
    @syms I_n real=true
    @syms omega domega theta_PLL dtheta_PLL v_xm dvx_m v_ym dv_ym P dP real=true
    @syms f_INVALID k_PLL_INVALID negative=true

    @test_throws AssertionError construct_vertex(PVInverterWithFrequencyControl(I_n=I_n,k_PLL=k_PLL_INVALID,f=f,f_s=f_s,T_m=T_m,k_P=k_P))
    @test_throws AssertionError construct_vertex(PVInverterWithFrequencyControl(I_n=I_n,k_PLL=k_PLL,f=f_INVALID,f_s=f_s,T_m=T_m,k_P=k_P))

    pv_inverter = PVInverterWithFrequencyControl(I_n=I_n,k_PLL=k_PLL,f=f,f_s=f_s,T_m=T_m,k_P=k_P)
    pv_inverter_vertex = construct_vertex(pv_inverter)

    @test symbolsof(pv_inverter) == [:u_r, :u_i,:θ_PLL,:v_xm,:v_ym,:P]
    @test pv_inverter_vertex.mass_matrix == [0,0,1,1,1,1]

    smoketest_rhs(pv_inverter_vertex, int_x=[omega,theta_PLL, v_xm, v_ym, P], int_dx=[domega, dtheta_PLL, dvx_m, dv_ym, dP])
end
