using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: SwingEq, construct_vertex, symbolsof
using LinearAlgebra: I

include("NodeTestBase.jl")

@testset "Swing Eq" begin
    @syms H D positive=true
    @syms P Ω real=true
    @syms omega domega real=true
    @syms H_INVALID D_INVALID negative=true

    @test_throws AssertionError construct_vertex(SwingEq(H=H, P=P, D=D_INVALID, Ω=Ω))
    @test_throws AssertionError construct_vertex(SwingEq(H=H_INVALID, P=P, D=D, Ω=Ω))

    swing = SwingEq(H=H, P=P, D=D, Ω=Ω)
    swing_vertex = construct_vertex(swing)

    @test symbolsof(swing) == [:u_r, :u_i, :ω]
    @test swing_vertex.mass_matrix == I

    smoketest_rhs(swing_vertex, int_x=omega, int_dx=domega)
end
