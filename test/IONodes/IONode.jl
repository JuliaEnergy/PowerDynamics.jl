using PowerDynamics
using BlockSystems
using PowerDynamics.IOComponents

@testset "BlockPara testes" begin
    blk = DroopControl(; x_ref=:a, K=:b, u_ref=:c)
    @testset "to_symbol" begin
        @parameters a t
        @variables x
        using PowerDynamics: to_symbol
        @test to_symbol(blk, blk.a) == :a
        @test to_symbol(blk, :a) == :a
        @test to_symbol(blk, a) == :a
        @test to_symbol(blk, BlockSystems.remove_namespace(blk.a)) == :a

        @test to_symbol(blk, blk.x) == :x
        @test to_symbol(blk, :x) == :x
        @test to_symbol(blk, x) == :x
        @test to_symbol(blk, BlockSystems.remove_namespace(blk.x)) == :x
    end

    @testset "Block para constructors" begin
        @parameters c
        p = Dict(blk.a=>1, :b=>2, c=>3)
        BlockPara(blk, p)

        p = Dict(blk.a=>1, :b=>2, c=>3, :d=>4, :e=>5)
        # additional p
        @test_throws ArgumentError BlockPara(blk, p)
        @test_throws ArgumentError BlockPara(blk, p; strict=false)

        p = Dict(blk.a=>1, :b=>2)
        @test_throws ArgumentError BlockPara(blk, p)
        @test_nowarn BlockPara(blk, p; strict=false)
    end

    @testset "mergep testes" begin
        using PowerDynamics: mergep

        @named blkA = DroopControl(; x_ref=:a, K=:b, u_ref=:c)
        pA = Dict(:a=>1, :b=>2, :c=>3)
        bp_A = BlockPara(blkA, pA)

        @named blkB = VoltageSource()
        pB = Dict(:τ=>4)
        bp_B = BlockPara(blkB, pB)

        # just one
        blks, paras = mergep(bp_A)
        @test blks == (blkA, )
        @test paras == Dict(:blkA₊a=>1, :blkA₊b=>2, :blkA₊c=>3)

        # merge two
        blks, paras = mergep(bp_A, bp_B)
        @test blks == (blkA, blkB)
        @test paras == Dict(:blkA₊a=>1, :blkA₊b=>2, :blkA₊c=>3, :blkB₊τ=>4)
    end
end

@testset "IONode constructor testes" begin
    blk = PowerConstraint()
    p = Dict(blk.P => 1.0, :Q => 2)
    node1 = IONode(blk, p);

    @test node1.parameter_names == [:P, :Q]
    @test node1.parameters == [1.0, 2.0]
    @test node1.mass_matrix == [0 0; 0 0]

    blk = DroopControl()
    p = Dict(:x_ref=>1.0, :K=>2.0, :u_ref=>3.0)
    # does not fulfil i/o convention
    @test_throws ArgumentError IONode(blk, p)
end
