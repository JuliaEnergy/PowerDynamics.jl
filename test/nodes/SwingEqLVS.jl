using Test
using SymPy
using PowerDynBase
using LinearAlgebra

include("NodeTestBase.jl")

@testset "Swing Eq LVS" begin
    @syms H D positive=true
    @syms P Ω real=true
    @syms omega domega real=true
    @syms V Γ positive=true
    swing_lvs = SwingEqLVS(H=H, P=P, D=D, Ω=Ω, Γ=Γ, V=V)
    swing_lvs_vertex = construct_vertex(swing_lvs)

    @test symbolsof(swing_lvs) == [:u_r, :u_i, :ω]
    @test swing_lvs_vertex.mass_matrix == I

    smoketest_rhs(swing_lvs_vertex, int_x=omega, int_dx=domega)
end
