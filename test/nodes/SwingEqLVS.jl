using Test: @testset, @test
using SymPy: @syms
using PowerDynBase: SwingEqLVS, construct_vertex, symbolsof
using LinearAlgebra: I

include("NodeTestBase.jl")

@testset "Swing Eq LVS" begin
    @syms H D positive=true
    @syms P Ω real=true
    @syms omega domega real=true
    @syms V Γ positive=true
    @syms H_INVALID D_INVALID Γ_INVALID negative=true
    swing_lvs = SwingEqLVS(H=H, P=P, D=D, Ω=Ω, Γ=Γ, V=V)
    swing_lvs_vertex = construct_vertex(swing_lvs)

    @test_throws AssertionError construct_vertex(SwingEqLVS(H=H_INVALID, P=P, D=D, Ω=Ω, Γ=Γ, V=V))
    @test_throws AssertionError construct_vertex(SwingEqLVS(H=H, P=P, D=D_INVALID, Ω=Ω, Γ=Γ, V=V))
    @test_throws AssertionError construct_vertex(SwingEqLVS(H=H, P=P, D=D, Ω=Ω, Γ=Γ_INVALID, V=V))

    @test symbolsof(swing_lvs) == [:u_r, :u_i, :ω]
    @test swing_lvs_vertex.mass_matrix == I

    smoketest_rhs(swing_lvs_vertex, int_x=omega, int_dx=domega)

    @test convert(SwingEq, swing_lvs) == SwingEq(H=H, P=P, D=D, Ω=Ω)
end
