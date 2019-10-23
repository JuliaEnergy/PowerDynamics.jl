using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: FESS, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "FESS Tests" begin
    @syms J k_PLL ω_m_ref f C L_g L_m P_n R_m R_g K_m1 K_m2 K_m3 K_m4 K_g1 K_g2 K_g3 K_g4 positive=true
    @syms u_dc_ref Ψ_m real=true
    @syms θ_PLL dθ_PLL i_dm di_dm i_qm di_qm ω_m dω_m v_qm dv_qm i_qm_ref di_qm_ref i_dg di_dg i_qg di_qg u_dc du_dc u_dg du_dg u_qg du_qg i_dg_ref di_dg_ref real=true
    @syms f_INVALID k_PLL_INVALID negative=true

    @test_throws AssertionError construct_vertex(FESS(u_dc_ref=u_dc_ref,C=C,J=J,ω_m_ref=ω_m_ref,k_PLL=k_PLL_INVALID,f=f,L_g=L_g,L_m=L_m,P_n=P_n,Ψ_m=Ψ_m,R_m=R_m,R_g=R_g,K_m1=K_m1,K_m2=K_m2,K_m3=K_m3,K_m4=K_m4,K_g1=K_g1,K_g2=K_g2,K_g3=K_g3,K_g4=K_g4))
    @test_throws AssertionError construct_vertex(FESS(u_dc_ref=u_dc_ref,C=C,J=J,ω_m_ref=ω_m_ref,k_PLL=k_PLL,f=f_INVALID,L_g=L_g,L_m=L_m,P_n=P_n,Ψ_m=Ψ_m,R_m=R_m,R_g=R_g,K_m1=K_m1,K_m2=K_m2,K_m3=K_m3,K_m4=K_m4,K_g1=K_g1,K_g2=K_g2,K_g3=K_g3,K_g4=K_g4))

    fess = FESS(u_dc_ref=u_dc_ref,C=C,J=J,ω_m_ref=ω_m_ref,k_PLL=k_PLL,f=f,L_g=L_g,L_m=L_m,P_n=P_n,Ψ_m=Ψ_m,R_m=R_m,R_g=R_g,K_m1=K_m1,K_m2=K_m2,K_m3=K_m3,K_m4=K_m4,K_g1=K_g1,K_g2=K_g2,K_g3=K_g3,K_g4=K_g4)
    fess_vertex = construct_vertex(fess)

    @test symbolsof(fess) == [:u_r, :u_i,:θ_PLL,:i_dm,:i_qm,:ω_m,:v_qm,:i_qm_ref,:i_dg,:i_qg,:u_dc,:u_dg,:u_qg,:i_dg_ref]
    @test fess_vertex.mass_matrix == [0,0,1,1,1,1,1,1,1,1,1,1,1,1]

    smoketest_rhs(fess_vertex, int_x=[θ_PLL,i_dm,i_qm,ω_m,v_qm,i_qm_ref,i_dg,i_qg,u_dc,u_dg,u_qg,i_dg_ref], int_dx=[dθ_PLL,di_dm,di_qm,dω_m,dv_qm,di_qm_ref,di_dg,di_qg,du_dc,du_dg,du_qg,di_dg_ref])
end
