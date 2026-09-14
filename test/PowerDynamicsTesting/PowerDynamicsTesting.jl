module PowerDynamicsTesting

# transitive dependencies (less restricitons on the env which loads this mod)
using PowerDynamics: OrderedDict
using PowerDynamics.NetworkDynamics.Graphs: path_graph
using PowerDynamics: Network, NWState, uflat, pflat
using PowerDynamics.NetworkDynamics: VIndex, EIndex, SII, VertexModel,
                                     ComponentAffect, PresetTimeComponentCallback, set_callback!
using PowerDynamics.SciMLBase: SciMLBase, solve, ODEProblem, auto_dt_reset!
using PowerDynamics.NetworkDynamics.DiffEqCallbacks: PresetTimeCallback
using PowerDynamics.ModelingToolkitBase: @named

using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve, BrownFullBasicInit
using Makie: Makie, Figure, Axis, axislegend, lines!, Cycled

using PowerDynamics: PowerDynamics, compile_bus, compile_line, MTKLine, initialize_from_pf!
using PowerDynamics.Library: PiLine, SlackDifferential, PSSE_Load, PSSE_GENCLS, PSSE_GENROE
using LinearAlgebra: norm
using Statistics: mean

using JLD2: JLD2
using Test: Test, @test, @test_broken

export TrajectoriesOfInterest, plottoi, compare
include("TrajectoriesOfInterest.jl")

export line_between_slacks, bus_on_slack
include("scenarios.jl")

export @reftest, set_reference_dir, refup
include("reftests.jl")

using PowerDynamics.Library: SauerPaiMachine, ConstantYLoad, AVRTypeI, TGOV1
using PowerDynamics: CompositeInjector, MTKBus, pfSlack, pfPV, pfPQ
using PowerDynamics: set_ωbase!, set_fbase!, get_ωbase
include("testsystems.jl")

export OpenIPSL_SMIB, ref_rms_error
include("OpenIPSLUtils.jl")

using IOCapture: IOCapture

"""
Include a test file into `mod` and throw its printed output away, unless something in the file
failed. Otherwise ParallelTestRunner echoes the output of every single file after the run.
"""
function quiet_include(mod, path)
    c = IOCapture.capture(; rethrow=Union{}, color=true) do
        Base.include(mod, path)
    end
    (c.error || _anynonpass(Test.get_testset())) && print(c.output)
    c.error && throw(CapturedException(c.value, c.backtrace))
    nothing
end
_anynonpass(res) = res isa Test.Fail || res isa Test.Error
_anynonpass(ts::Test.AbstractTestSet) = any(_anynonpass, ts.results)

end # module PowerDynamicsTesting
