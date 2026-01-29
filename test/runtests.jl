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
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

@testset "PowerDynamics.jl Tests" begin
    @testset "Package Quality Tests" begin
        @info "Begin Package quatlity tests"
        Aqua.test_all(PowerDynamics;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(PowerDynamics))

        MTKMODEL_SYMS = (:Num, :System, :Equation, :connect, :@unpack, :setmetadata, :ComponentPostprocessing)
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

    @safetestset "Library tests" begin include("Library_test.jl") end
    @safetestset "Saturation tests" begin include("saturation_test.jl") end
    @safetestset "utils tests" begin include("utils_test.jl") end
    @safetestset "modeling_tools tests" begin include("modeling_tools_test.jl") end
    @safetestset "initialization tests" begin include("initialization_test.jl") end

    EXPORT_FIGURES = false
    @testset "OpenIPSL Model Tests" begin
        # Machines
        @safetestset "PSSE_GENCLS" begin include(joinpath("OpenIPSL_test", "PSSE_GENCLS_test.jl")) end
        @safetestset "PSSE_GENROU" begin include(joinpath("OpenIPSL_test", "PSSE_GENROU_test.jl")) end
        @safetestset "PSSE_GENROE" begin include(joinpath("OpenIPSL_test", "PSSE_GENROE_test.jl")) end
        @safetestset "PSSE_GENSAL" begin include(joinpath("OpenIPSL_test", "PSSE_GENSAL_test.jl")) end
        @safetestset "PSSE_GENSAE" begin include(joinpath("OpenIPSL_test", "PSSE_GENSAE_test.jl")) end

        # Exciters
        @safetestset "PSSE_IEEET1" begin include(joinpath("OpenIPSL_test", "PSSE_IEEET1_test.jl")) end
        @safetestset "PSSE_SCRX" begin include(joinpath("OpenIPSL_test", "PSSE_SCRX_test.jl")) end
        @safetestset "PSSE_ESST4B" begin include(joinpath("OpenIPSL_test", "PSSE_ESST4B_test.jl")) end
        @safetestset "PSSE_EXST1" begin include(joinpath("OpenIPSL_test", "PSSE_EXST1_test.jl")) end
        @safetestset "PSSE_ESST1A" begin include(joinpath("OpenIPSL_test", "PSSE_ESST1A_test.jl")) end

        # Governors
        @safetestset "PSSE_IEEEG1" begin include(joinpath("OpenIPSL_test", "PSSE_IEEEG1_test.jl")) end
        @safetestset "PSSE_GGOV1" begin include(joinpath("OpenIPSL_test", "PSSE_GGOV1_test.jl")) end
        @safetestset "PSSE_HYGOV" begin include(joinpath("OpenIPSL_test", "PSSE_HYGOV_test.jl")) end

        # Power System Stabilizers
        @safetestset "PSSE_IEEEST" begin include(joinpath("OpenIPSL_test", "PSSE_IEEEST_test.jl")) end
    end

    @testset "validation tests" begin
        @safetestset "ieee39 RMSPowerSims.jl" begin include(joinpath("validation", "ieee39_RMSPowerSims.jl", "ieee39_validation.jl")) end
    end

    @testset "Test Doc Tutorials" begin
        examples = joinpath(pkgdir(PowerDynamics), "docs", "tutorials")
        for file in readdir(examples; join=true)
            endswith(file, ".jl") || continue
            name = basename(file)
            @info "Test Tutorial $name"
            eval(:(@safetestset $name begin include($file) end))
        end
    end

    @testset "Test Doc Examples" begin
        examples = joinpath(pkgdir(PowerDynamics), "docs", "examples")
        for file in readdir(examples; join=true)
            endswith(file, ".jl") || continue
            name = basename(file)
            @info "Test Example $name"
            eval(:(@safetestset $name begin include($file) end))
        end
    end
end
