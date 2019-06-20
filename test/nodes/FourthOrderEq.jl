using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: FourthOrderEq, construct_vertex, symbolsof
using LinearAlgebra: I

include("NodeTestBase.jl")

@testset "FourthOrderEq" begin
    @syms H  D  Ω  T_d_dash T_q_dash X_q_dash X_d_dash X_d X_q positive=true
    @syms P  E_f real=true
    @syms omega domega theta dtheta real=true
    fourth_order = FourthOrderEq(H=H, P=P, D=D, Ω=Ω, E_f=E_f, T_d_dash=T_d_dash ,T_q_dash=T_q_dash ,X_q_dash=X_q_dash ,X_d_dash=X_d_dash,X_d=X_d, X_q=X_q)
    fourth_order_vertex = construct_vertex(fourth_order)

    @test symbolsof(fourth_order) == [:u_r, :u_i, :θ, :ω]
    @test fourth_order_vertex.mass_matrix == I

    smoketest_rhs(fourth_order_vertex, int_x=[theta, omega], int_dx=[dtheta, domega])
end
