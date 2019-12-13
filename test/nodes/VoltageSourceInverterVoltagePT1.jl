using Test: @testset, @test
using PowerDynamics: VSIVoltagePT1, construct_vertex, symbolsof
using LinearAlgebra: I

@testset "VSIVoltagePT1" begin
    τ_v, τ_P, τ_Q, K_P, K_Q = rand_positive(5)
    P,Q,V_r = rand_real(3)
    q_m, dq_m, omega, domega = rand_real(4)
    VSI_PT1 = VSIVoltagePT1(τ_v=τ_v,τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)
    VSI_PT1_vertex = construct_vertex(VSI_PT1)
    dint = [domega,dq_m]; int = [omega,q_m]; int_test = copy(int)


    @test symbolsof(VSI_PT1) == [:u_r, :u_i, :ω, :q_m]
    @test VSI_PT1_vertex.mass_matrix == I

    smoketest_rhs(VSI_PT1_vertex, int_x=[q_m, omega], int_dx=[dq_m, domega])
end
