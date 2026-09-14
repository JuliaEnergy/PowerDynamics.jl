using PowerDynamics
using PowerDynamics.Library
using NetworkDynamics
using Main.PowerDynamicsTesting
using Aqua
using ExplicitImports
using Test

@testset "Package Quality Tests" begin
    Aqua.test_all(PowerDynamics;
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
