using Test: @testset, @test
using PowerDynamics: FourthOrderEqGovernorIEEEG1, construct_vertex, symbolsof
using LinearAlgebra: I

include("NodeTestBase.jl")

@testset "FourthOrderEqGovernorIEEEG1" begin
    H, D, Ω, E_f, T_d_dash, T_q_dash, X_d_dash, X_q_dash, X_d, X_q, P0, Pmax, Pmin, Pup, Pdown, T1, T2, T3, K = rand_positive(19)
    omega, domega, theta, dtheta, Pm, dPm, x1, dx1, z, dz, P, dP = rand_real(12)
    fourth_order = FourthOrderEqGovernorIEEEG1(H = H, D = D, Ω = Ω, E_f = E_f, 
                                               T_d_dash = T_d_dash,  T_q_dash = T_q_dash, X_d_dash = X_d_dash,  X_q_dash = X_q_dash, X_d = X_d, X_q = X_q, 
                                               P0 = P0, Pmax = Pmax, Pmin = Pmin, Pup = Pup, Pdown = Pdown, T1 = T1, T2 = T2, T3 = T3, K = K)
    fourth_order_vertex = construct_vertex(fourth_order)

    @test symbolsof(fourth_order) == [:u_r, :u_i, :θ, :ω, :Pm, :x1, :z, :P]
    @test fourth_order_vertex.mass_matrix == I

    smoketest_rhs(fourth_order_vertex, int_x = [theta, omega, Pm, x1, z, P], int_dx = [dtheta, domega, dPm, dx1, dz, dP])
end
