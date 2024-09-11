using Test
using SafeTestsets
using Aqua
using ExplicitImports
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using Makie
# using GLMakie
using OrderedCollections

using OpPoDyn
using OpPoDyn.Library
using OpPoDynTesting
set_reference_dir(OpPoDyn)

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


    @testset "pin parameters" begin
        @mtkmodel Inner begin
            @parameters begin
                a = 0
                b = 0
            end
            @variables begin
                i(t)
            end
            @equations begin
                i ~ a + b
            end
        end
        @mtkmodel Outer begin
            @components begin
                inner = Inner()
            end
            @parameters begin
                a = 0
                b = 0
            end
            @variables begin
                o(t)
            end
            @equations begin
                o ~ inner.i * (a + b)
            end
        end
        @named outer = Outer()
        full_equations(outer)
        sys1 = pin_parameters(outer, :a => 1)
        sys2 = pin_parameters(outer, :outer₊b => 1)
        full_equations(sys1)
        full_equations(sys2)
        sys1 = pin_parameters(outer, outer.inner.a => 1)
        sys2 = pin_parameters(outer, :inner₊a => 1)
        @test full_equations(sys1) == full_equations(sys2)
    end
end
