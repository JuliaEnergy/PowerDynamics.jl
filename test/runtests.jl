using Test
using SafeTestsets
using Aqua
using ExplicitImports
using NetworkDynamics
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using Makie
using CairoMakie
using OrderedCollections

using OpPoDyn
using OpPoDyn.Library
using OpPoDynTesting

@testset "OpPoDyn.jl Tests" begin
    @testset "Package Quality Tests" begin
        Aqua.test_all(OpPoDyn;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(OpPoDyn))

        @test check_no_implicit_imports(OpPoDyn) === nothing
        @test_broken check_no_stale_explicit_imports(OpPoDyn) === nothing

        path = joinpath(pkgdir(OpPoDyn),"src","Library","Library.jl")
        @test check_no_implicit_imports(OpPoDyn.Library, path) === nothing
        @test_broken check_no_stale_explicit_imports(OpPoDyn.Library, path) === nothing


        # check the testing subpackage
        Aqua.test_all(OpPoDynTesting;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(OpPoDynTesting))

        @test check_no_implicit_imports(OpPoDynTesting) === nothing
        @test_broken check_no_stale_explicit_imports(OpPoDynTesting) === nothing
    end

    @safetestset "Library tests" begin include("Library_test.jl") end
    @safetestset "utils tests" begin include("utils_test.jl") end
    @safetestset "initialization tests" begin include("initialization_test.jl") end

    @testset "validation tests" begin
        @safetestset "ieee39 RMSPowerSims.jl" begin include(joinpath("validation", "ieee39_RMSPowerSims.jl", "ieee39_validation.jl")) end
    end
end
