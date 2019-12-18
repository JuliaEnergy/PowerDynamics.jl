using Test: @testset, @test
using PowerDynamics: FourthOrderEqGovernorExciterAVR, construct_vertex, symbolsof
using LinearAlgebra: I

@testset "FourthOrderEqGovernorExciterAVR" begin
    H,D,Ω,T_d_dash,T_q_dash,X_q_dash,X_d_dash,X_d,X_q,T_e,T_a,T_f,K_a,K_f,V_ref,R_d,T_sv,T_ch = rand_positive(18)
    P,K_e = rand_real(2)
    omega,domega,theta,dtheta,e_f,de_f,v_r,dv_r,r_f,dr_f,P_sv,dP_sv,P_m,dP_m = rand_real(14)
    fourth = FourthOrderEqGovernorExciterAVR(
    H=H, P=P, D=D, Ω=Ω, T_d_dash=T_d_dash ,T_q_dash=T_q_dash,
    X_q_dash=X_q_dash ,X_d_dash=X_d_dash,X_d=X_d, X_q=X_q, T_e=T_e, T_a=T_a, T_f=T_f, K_e=K_e, K_a=K_a,
    K_f=K_f, V_ref=V_ref, R_d=R_d, T_sv=T_sv, T_ch=T_ch)

    fourth_vertex = construct_vertex(fourth)

    @test symbolsof(fourth) == [:u_r, :u_i, :θ, :ω, :e_f, :v_r, :r_f, :P_sv, :P_m]
    @test fourth_vertex.mass_matrix == I

    smoketest_rhs(fourth_vertex,
    int_x=[theta,omega,e_f,v_r,r_f,P_sv,P_m],
    int_dx=[dtheta,domega,de_f,dv_r,dr_f,dP_sv,dP_m])
end
