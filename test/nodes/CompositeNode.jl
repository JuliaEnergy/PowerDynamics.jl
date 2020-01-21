using Revise

using Test: @testset, @test
using PowerDynamics#: CompositeNode, CSIMinimal, VSIMinimal, PQAlgebraic, construct_vertex, symbolsof

# constant voltage
VSI = VSIMinimal(τ_P=1.,τ_Q=1.,K_P=1.,K_Q=1.,V_r=1.,P=1.,Q=1.)

# constant current
CSI = CSIMinimal(I_r=1.)
CSI2 = CSIMinimal(I_r=1im)

# constant power
PQ = PQAlgebraic(P=0.1,Q=0.1)
PQ2 = PQAlgebraic(P=0.2,Q=0.2)

comp_node = CompositeNode(CurrentNodes=[CSI, CSI2], PowerNodes=[PQ, PQ2], VoltageNode=VSI)

typeof(comp_node)

construct_vertex(comp_node)


@testset "CompositeNode" begin
    # constant voltage
    VSI = VSIMinimal(τ_P=1.,τ_Q=1.,K_P=1.,K_Q=1.,V_r=1.,P=1.,Q=1.)

    # constant current
    CSI = CSIMinimal(I_r=1.)
    CSI2 = CSIMinimal(I_r=2.)

    # constant power
    PQ = PQAlgebraic(P=0.1,Q=0.1)
    PQ2 = PQAlgebraic(P=0.2,Q=0.2)

    comp_node = CompositeNode(CurrentNodes=[CSI, CSI2], PowerNodes=[PQ, PQ2], VoltageNode=VSI)
    comp_vertex = construct_vertex(comp_node)

    @test symbolsof(comp_node) == [:u_r, :u_i]
    @test comp_vertex.mass_matrix == spzeros(2,2)

    smoketest_rhs(comp_vertex)
end
