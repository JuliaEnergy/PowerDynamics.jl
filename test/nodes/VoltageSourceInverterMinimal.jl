using Test: @testset, @test
using PowerDynamics: VSIMinimal, construct_vertex, symbolsof
using LinearAlgebra: I

@testset "VoltageSourceInverterMinimal" begin
    τ_P,τ_Q,K_P,K_Q = rand_positive(4)
    P,Q,V_r = rand_real(3)
    omega,domega = rand_real(2)
    VSI = VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)
    VSI_vertex = construct_vertex(VSI)

    @test symbolsof(VSI) == [:u_r, :u_i, :ω]
    @test VSI_vertex.mass_matrix == I

    smoketest_rhs(VSI_vertex, int_x=[omega], int_dx=[domega])
end
