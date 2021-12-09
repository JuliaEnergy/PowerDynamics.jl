using Test: @testset, @test
using PowerDynamics.IOComponents
using PowerDynamics: BlockPara, BusNode, symbolsof, construct_vertex
using LinearAlgebra: Diagonal

@testset "BusNode construction" begin
    pq_constraint = InversePowerConstraint(;name=:inv_pqconstraint,P=:P,Q=:Q)
    pq_para = Dict(pq_constraint.P => rand(),
                pq_constraint.Q => rand())
    pq_load = BlockPara(pq_constraint, pq_para)

    rx_constraint = ImpedanceConstraint(;name=:rxconstraint)
    rx_para = Dict(rx_constraint.R => rand(),
                rx_constraint.X => rand())
    rx_load = BlockPara(rx_constraint, rx_para)

    busnode = BusNode(pq_load,rx_load)

    @test symbolsof(busnode) == [:u_r, :u_i]
    @test busnode.mass_matrix == Diagonal([0,0])

    busnode_vertex = construct_vertex(busnode)
    smoketest_rhs(busnode_vertex)

    emptybus = BusNode()

    @test symbolsof(emptybus) == [:u_r, :u_i]
    @test emptybus.mass_matrix == Diagonal([0,0])

    emptybus_vertex = construct_vertex(emptybus)
    smoketest_rhs(emptybus_vertex)
end