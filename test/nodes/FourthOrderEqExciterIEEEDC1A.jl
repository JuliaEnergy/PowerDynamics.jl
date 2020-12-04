using Test: @testset, @test
using PowerDynamics: FourthOrderEqExciterIEEED1A, construct_vertex, symbolsof
using LinearAlgebra: I

include("NodeTestBase.jl")

@testset "FourthOrderEqExciterIEEED1A" begin
    H, P, D, Ω, T_d_dash, T_q_dash, X_d_dash, X_q_dash, X_d, X_q, K_e, K_f, K_a, U, U_ref, U_ref2, U_rmax, U_rmin, T_a, T_f, T_e, S_E_max, S_E_tq, V_R_max = rand_positive(24)
   
    theta, dtheta, omega, domega, E_f, dE_f, U_f, dU_f, U_r, dU_r = rand_real(10)
    fourth_order = FourthOrderEqExciterIEEED1A(H=H, P=P, D=D, Ω=Ω, 
                                               T_d_dash=T_d_dash, T_q_dash=T_q_dash, X_d_dash=X_d_dash, X_q_dash=X_q_dash, X_d=X_d, X_q=X_q, 
                                               K_e=K_e, K_f=K_f, K_a=K_a, U=U, U_ref=U_ref, U_ref2=U_ref2, U_rmax=U_rmax, U_rmin=U_rmin, T_a=T_a, T_f=T_f, T_e=T_e, S_E_max=S_E_max, S_E_tq=S_E_max, V_R_max=V_R_max)
    fourth_order_vertex = construct_vertex(fourth_order)

    @test symbolsof(fourth_order) == [:u_r, :u_i, :θ, :ω, :E_f, :U_f, :U_r]
    @test fourth_order_vertex.mass_matrix == I

    smoketest_rhs(fourth_order_vertex, int_x=[theta, omega, E_f, U_f, U_r], int_dx=[dtheta, domega, dE_f, dU_f, dU_r])
end
