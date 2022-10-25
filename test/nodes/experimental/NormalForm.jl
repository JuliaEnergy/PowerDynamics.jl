using Test: @testset, @test
using PowerDynamics: NormalForm, construct_vertex, symbolsof, dimension
using LinearAlgebra: diag, I

include("../NodeTestBase.jl")

@testset "NormalForm Tests" begin

    P, Q, V = rand_real(3)

    # Case 1: no internal variable
    Cᵤ, Gᵤ, Hᵤ = 1im*rand_real(3)

    nf = NormalForm(P=P, Q=Q, V=V, Cᵤ=Cᵤ, Gᵤ=Gᵤ, Hᵤ=Hᵤ)
    nf_vertex = construct_vertex(nf)

    @test dimension(nf) == 2
    @test symbolsof(nf) == [:u_r, :u_i]
    @test nf_vertex.mass_matrix == I

    smoketest_rhs(nf_vertex)

    # Case 2: one internal variable
    Bᵤ, Cᵤ, Gᵤ, Hᵤ = 1im*rand_real(4)
    Bₓ, Cₓ, Gₓ, Hₓ = rand_real(4)

    nf = NormalForm(P=P, Q=Q, V=V, Bᵤ=Bᵤ, Cᵤ=Cᵤ, Gᵤ=Gᵤ, Hᵤ=Hᵤ, Bₓ=Bₓ, Cₓ=Cₓ, Gₓ=Gₓ, Hₓ=Hₓ)
    nf_vertex = construct_vertex(nf)

    @test dimension(nf) == 3
    @test symbolsof(nf) == [:u_r, :u_i, :x_1]
    @test nf_vertex.mass_matrix == I

    x_1, dx_1 = rand_real(2)
    smoketest_rhs(nf_vertex, int_x=[x_1], int_dx=[dx_1])

    # Case 3: two internal variables
    Bᵤ = 1im*rand_real(2)
    Cᵤ, Gᵤ, Hᵤ = 1im*rand_real(3)
    Bₓ = rand(2,2)
    Cₓ = rand_real(2)
    Gₓ = rand_real(2)
    Hₓ = rand_real(2)

    nf = NormalForm(P=P, Q=Q, V=V, Bᵤ=Bᵤ, Cᵤ=Cᵤ, Gᵤ=Gᵤ, Hᵤ=Hᵤ, Bₓ=Bₓ, Cₓ=Cₓ, Gₓ=Gₓ, Hₓ=Hₓ)
    nf_vertex = construct_vertex(nf)

    @test dimension(nf) == 4
    @test symbolsof(nf) == [:u_r, :u_i, :x_1, :x_2]
    @test nf_vertex.mass_matrix == I

    x_1, dx_1, x_2, dx_2 = rand_real(4)
    smoketest_rhs(nf_vertex, int_x=[x_1,x_2], int_dx=[dx_1,dx_2])

end