using ModelingToolkit
using Test
using Testfiles
using Aqua
using ExplicitImports
using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using Makie
using CairoMakie
using OrderedCollections

using PowerDynamics
using PowerDynamics.Library
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

@testset "PowerDynamics.jl Tests" begin
    @testset "Package Quality Tests" begin
        @info "Begin Package quatlity tests"
        Aqua.test_all(PowerDynamics;
            deps_compat=VERSION ≥ v"1.11", # don't check compat on LTS
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(PowerDynamics))

        MTKMODEL_SYMS = (:Num, :System, :Equation, :connect, Symbol("@unpack"), :setmetadata, :ComponentPostprocessing, Symbol("@named"))
        allow_unanalyzable = (PowerDynamics,)

        @test check_no_implicit_imports(PowerDynamics; skip=(Base, Core, NetworkDynamics), allow_unanalyzable) === nothing
        @test check_no_stale_explicit_imports(PowerDynamics; ignore=MTKMODEL_SYMS, allow_unanalyzable) === nothing

        path = joinpath(pkgdir(PowerDynamics),"src","Library","Library.jl")
        @test check_no_implicit_imports(PowerDynamics.Library, path) === nothing
        @test check_no_stale_explicit_imports(PowerDynamics.Library, path; ignore=MTKMODEL_SYMS) === nothing

        pdt_path = PowerDynamics.pdtesting_path()
        @test check_no_implicit_imports(PowerDynamicsTesting, pdt_path) === nothing
        @test check_no_stale_explicit_imports(PowerDynamicsTesting, pdt_path) === nothing
    end

    @testfile "Library_test.jl"
    @testfile "saturation_test.jl"
    @testfile "utils_test.jl"
    @testfile "modeling_tools_test.jl"
    @testfile "initialization_test.jl"

    EXPORT_FIGURES = false
    @testset "OpenIPSL Model Tests" begin
        # Machines
        @testfile joinpath("OpenIPSL_test", "PSSE_GENCLS_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_GENROU_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_GENROE_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_GENSAL_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_GENSAE_test.jl")

        # Exciters
        @testfile joinpath("OpenIPSL_test", "PSSE_IEEET1_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_SCRX_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_ESST4B_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_EXST1_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_ESST1A_test.jl")

        # Governors
        @testfile joinpath("OpenIPSL_test", "PSSE_IEEEG1_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_GGOV1_test.jl")
        @testfile joinpath("OpenIPSL_test", "PSSE_HYGOV_test.jl")

        # Power System Stabilizers
        @testfile joinpath("OpenIPSL_test", "PSSE_IEEEST_test.jl")
    end

    @testset "validation tests" begin
        @testfile joinpath("validation", "ieee39_RMSPowerSims.jl", "ieee39_validation.jl")
    end

    @testset "Test Doc Tutorials" begin
        examples = joinpath(pkgdir(PowerDynamics), "docs", "tutorials")
        for file in readdir(examples; join=true)
            endswith(file, ".jl") || continue
            eval(:(@testfile $file))
        end
    end

    @testset "Test Doc Examples" begin
        examples = joinpath(pkgdir(PowerDynamics), "docs", "examples")
        for file in readdir(examples; join=true)
            endswith(file, ".jl") || continue
            eval(:(@testfile $file))
        end
    end
end
