using Test: @testset, @test
using SymPy: @syms
using PowerDynBase: CSIMinimal, construct_vertex, symbolsof

@testset "CSIMinimal" begin
    @syms I_r real=true
    CSI_min = CSIMinimal(I_r=I_r)
    CSI_min_vertex= construct_vertex(CSI_min)

    @test symbolsof(CSI_min) == [:u_r, :u_i]
    @test CSI_min_vertex.mass_matrix == [0,0]

    smoketest_rhs(CSI_min_vertex)
end
