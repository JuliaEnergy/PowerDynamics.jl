using Test: @testset, @test
using SymPy: @syms
using PowerDynBase: VSIMinimal, construct_vertex, symbolsof
using LinearAlgebra: I

@testset "VoltageSourceInverterMinimal" begin
    @syms τ_P τ_Q K_P K_Q positive=true
    @syms P Q V_r real=true
    @syms omega domega real=true
    VSI = VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)
    VSI_vertex = construct_vertex(VSI)

    @test symbolsof(VSI) == [:u_r, :u_i, :ω]
    @test VSI_vertex.mass_matrix == I

    smoketest_rhs(VSI_vertex, int_x=[omega], int_dx=[domega])
end
