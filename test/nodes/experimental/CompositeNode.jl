using Test: @testset, @test
using PowerDynamics: CompositeNode, CSIMinimal, VSIMinimal, PQAlgebraic, construct_vertex, symbolsof
using LinearAlgebra: diag

include("../NodeTestBase.jl")

@testset "CompositeNode" begin
    # constant voltage
    #VSI = VSIMinimal(τ_P=1.,τ_Q=1.,K_P=1.,K_Q=1.,V_r=1.,P=1.,Q=1.)

    τ_P, τ_Q, K_P, K_Q = rand_real(4)
    P,Q,V_r =rand_real(3)
    omega,domega = rand_real(2)
    VSI = VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)

    # constant current
    I_r =rand(1.0:10.0)+1im*rand(1.0:10.0)
    CSI = CSIMinimal(I_r=I_r)
    CSI2 = CSIMinimal(I_r=I_r)

    # constant power
    P,Q=rand_real(2)
    PQ = PQAlgebraic(P=P,Q=Q)
    PQ2 = PQAlgebraic(P=P,Q=Q)

    # Test vertex with voltage dynamic
    comp_node = CompositeNode(CurrentNodes=[CSI, CSI2], PowerNodes=[PQ, PQ2], VoltageNode=VSI)
    comp_vertex = construct_vertex(comp_node)

    @test symbolsof(comp_node) == [:u_r, :u_i, :ω]
    @test comp_vertex.mass_matrix |> diag == [1,1,1]

    smoketest_rhs(comp_vertex, int_x=omega, int_dx=domega)
    # smoketest not compatible with internal variables?
    # comp_vertex.f(ones(3), ones(3), [], [], nothing, 0.)

    # Test vertex without voltage dynamic
    comp_node_2 = CompositeNode(CurrentNodes=[CSI, CSI2], PowerNodes=[PQ, PQ2])
    comp_vertex_2 = construct_vertex(comp_node_2)

    @test symbolsof(comp_node_2) == [:u_r, :u_i]
    @test comp_vertex_2.mass_matrix |> diag == [0,0]

    smoketest_rhs(comp_vertex_2)
end
