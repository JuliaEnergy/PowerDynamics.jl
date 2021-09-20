using Test
using SafeTestsets
using PowerDynamicsNREL
using BlockSystems

@testset "BlockPara" begin
    using PowerDynamicsNREL: BlockPara, mergep

    b = PowerDynamicsNREL.base_machine()
    p1 = Dict(b.R => 1.0, b.X_d => 1.0, b.e_q => 1.0)
    p2 = Dict(b.R => 1.0, b.X_d => 1.0)
    p3 = Dict(b.R => 1.0, b.X_d => 1.0, b.e_q => 1.0, :foo => 2.0)

    @test BlockPara(b, p1) isa BlockPara
    @test_throws ArgumentError BlockPara(b, p2)
    @test BlockPara(b, p3) isa BlockPara

    # test merging
    b1 = PowerDynamicsNREL.base_machine()
    p1 = Dict(b1.R => 1.0, b1.X_d => 1.0, b1.e_q => 1.0)
    bp1 = BlockPara(b1, p1)

    b2 = PowerDynamicsNREL.avr_simple()
    p2 = Dict(b2.K => 2.0, b2.v_ref => 3.14)
    bp2 = BlockPara(b2, p2)

    b3 = PowerDynamicsNREL.single_mass_shaft()
    p3 = Dict(b3.Ω_b => 1.2, b3.ω_ref => 50, b3.D => 4, b3.H => 2.2)
    bp3 = BlockPara(b3, p3)

    @test_throws ArgumentError mergep(bp1, bp1)
    @test mergep(bp1, bp1; strict=false) isa Tuple{NTuple{2, IOBlock}, <:Dict}
    @test mergep(bp1, bp2, bp3) isa Tuple
end

@safetestset "MetaGenerator Tests" begin include("MetaGenerator_test.jl") end
@safetestset "Bus Tests" begin include("Bus_test.jl") end
