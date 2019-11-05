using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: WindTurbineGenType4, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "WindGenType 4 Tests" begin
    @syms ΔP_max K_P K_PLL positive=true
    @syms P real=true
    @syms θ_PLL dθ_PLL e_Idθ de_Idθ i_d di_d ω dω real=true
    @syms k_P_invalid k_PLL_INVALID negative=true

    @test_throws AssertionError construct_vertex(WindTurbineGenType4(ΔP_max=ΔP_max,K_P=K_P,K_PLL=k_PLL_INVALID,P=P))
    @test_throws AssertionError construct_vertex(WindTurbineGenType4(ΔP_max=ΔP_max,K_P=k_P_invalid,K_PLL=K_PLL,P=P))

    wind_gen = WindTurbineGenType4(ΔP_max=ΔP_max,K_P=K_P,K_PLL=K_PLL,P=P)
    wind_gen_vertex = construct_vertex(wind_gen)

    @test symbolsof(wind_gen) == [:u_r,:u_i,:θ_PLL,:ω,:e_Idθ,:i_d]
    @test wind_gen_vertex.mass_matrix == [0,0,1,1,1,1]

    smoketest_rhs(wind_gen_vertex, int_x=[θ_PLL,e_Idθ,i_d, ω], int_dx=[dθ_PLL,de_Idθ,di_d, dω])
end
