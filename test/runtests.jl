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

using PowerDynamics
using PowerDynamics.Library
using PowerDynamicsTesting

@testset "PowerDynamics.jl Tests" begin
    @testset "Package Quality Tests" begin
        @info "Begin Package quatlity tests"
        Aqua.test_all(PowerDynamics;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(PowerDynamics))

        @test check_no_implicit_imports(PowerDynamics; skip=(Base, Core, NetworkDynamics)) === nothing
        # mtkmodel macro depends on some symbols
        @test_broken check_no_stale_explicit_imports(PowerDynamics) === nothing

        path = joinpath(pkgdir(PowerDynamics),"src","Library","Library.jl")
        @test check_no_implicit_imports(PowerDynamics.Library, path) === nothing
        # mtkmodel macro depends on some symbols
        @test_broken check_no_stale_explicit_imports(PowerDynamics.Library, path) === nothing

        # check the testing subpackage
        Aqua.test_all(PowerDynamicsTesting;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(PowerDynamicsTesting))

        @test check_no_implicit_imports(PowerDynamicsTesting) === nothing
        @test_broken check_no_stale_explicit_imports(PowerDynamicsTesting) === nothing
    end

    @safetestset "Library tests" begin include("Library_test.jl") end
    @safetestset "utils tests" begin include("utils_test.jl") end
    @safetestset "modeling_tools tests" begin include("modeling_tools_test.jl") end
    @safetestset "initialization tests" begin include("initialization_test.jl") end

    @testset "validation tests" begin
        @safetestset "ieee39 RMSPowerSims.jl" begin include(joinpath("validation", "ieee39_RMSPowerSims.jl", "ieee39_validation.jl")) end
    end
end
