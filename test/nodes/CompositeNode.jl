using Test: @testset, @test
using PowerDynamics: CompositeNode, CSIMinimal, construct_vertex, symbolsof
using SparseArrays

@testset "CompositeNode" begin
    comp_node = CompositeNode([CSIMinimal(I_r = 1.0), CSIMinimal(I_r = 1.0im)])
    comp_vertex = construct_vertex(comp_node)

    @test symbolsof(comp_node) == [:u_r, :u_i]
    @test comp_vertex.mass_matrix == spzeros(2,2)

    smoketest_rhs(comp_vertex)
end
