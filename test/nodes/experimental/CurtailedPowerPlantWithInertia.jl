using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: CurtailedPowerPlantWithInertia, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "CurtailedPowerPlantWithInertia Tests" begin
    #T_L,T_H,K_P,K_PLL,Q_ref,C,J,P,ω_rref,u_dcref,K_Q,K_v,K_g1,K_g2,K_r1,K_r2
    #[θ_PLL,dθ_PLL],[e_Idθ,de_Idθ],[ω,dω],[e_IP,de_IP],[z,dz],[e_IV,de_IV],[u_dc,du_dc],[i_q,di_q],[u_tref,du_tref],[ω_r,dω_r]
    @syms ω_0 T_AI K_PPLL K_IPLL T_d T_f K_PV K_IV positive=true
    @syms P real=true
    @syms θ_PLL dθ_PLL e_Iω de_Iω ω dω y dy e_Iiq de_Iiq e_Iid de_Iid real=true
    @syms K_PPLL_invalid K_PV_INVALID negative=true

    @test_throws AssertionError construct_vertex(CurtailedPowerPlantWithInertia(P=P,ω_0=ω_0,T_AI=T_AI,K_PPLL=K_PPLL_invalid,K_IPLL=K_IPLL,T_d=T_d,T_f=T_f,K_PV=K_PV,K_IV=K_IV))
    @test_throws AssertionError construct_vertex(CurtailedPowerPlantWithInertia(P=P,ω_0=ω_0,T_AI=T_AI,K_PPLL=K_PPLL,K_IPLL=K_IPLL,T_d=T_d,T_f=T_f,K_PV=K_PV_INVALID,K_IV=K_IV))

    curtailed_pp = CurtailedPowerPlantWithInertia(P=P,ω_0=ω_0,T_AI=T_AI,K_PPLL=K_PPLL,K_IPLL=K_IPLL,T_d=T_d,T_f=T_f,K_PV=K_PV,K_IV=K_IV)
    curtailed_pp_vertex = construct_vertex(curtailed_pp)

    @test symbolsof(curtailed_pp) == [:u_r,:u_i,:θ_PLL,:e_Iω, :ω,:y,:e_Iiq,:e_Iid]
    @test curtailed_pp_vertex.mass_matrix == [0,0,1,1,1,1,1,1]

    smoketest_rhs(curtailed_pp_vertex, int_x=[θ_PLL,e_Iω, ω,y,e_Iiq,e_Iid], int_dx=[dθ_PLL,de_Iω, dω,dy,de_Iiq,de_Iid])
end
