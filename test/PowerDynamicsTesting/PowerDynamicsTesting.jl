module PowerDynamicsTesting

# transitive dependencies (less restricitons on the env which loads this mod)
using PowerDynamics: OrderedDict
using PowerDynamics.NetworkDynamics.Graphs: path_graph
using PowerDynamics: Network, NWState, uflat, pflat
using PowerDynamics.NetworkDynamics: VIndex, EIndex, SII, VertexModel,
                                     ComponentAffect, PresetTimeComponentCallback, set_callback!,
                                     get_callbacks
using PowerDynamics.SciMLBase: SciMLBase, solve, ODEProblem, auto_dt_reset!
using PowerDynamics.NetworkDynamics.DiffEqCallbacks: PresetTimeCallback
using PowerDynamics.ModelingToolkit: @named

using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
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
include("testsystems.jl")

export OpenIPSL_SMIB, ref_rms_error
include("OpenIPSLUtils.jl")

end # module PowerDynamicsTesting
