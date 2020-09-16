using Test: @testset, @test
using PowerDynamics: CSIMinimal, construct_vertex, symbolsof
using LinearAlgebra: diag

@testset "CSIMinimal" begin
    I_r = rand(1.0:10.0)
    CSI_min = CSIMinimal(I_r=I_r)
    CSI_min_vertex= construct_vertex(CSI_min)

    @test symbolsof(CSI_min) == [:u_r, :u_i]
    @test CSI_min_vertex.mass_matrix |> diag == [0,0]

    smoketest_rhs(CSI_min_vertex)
end
