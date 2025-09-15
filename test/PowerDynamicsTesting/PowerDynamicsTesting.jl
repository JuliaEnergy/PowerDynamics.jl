module PowerDynamicsTesting

# transitive dependencies (less restricitons on the env which loads this mod)
using PowerDynamics: OrderedDict
using PowerDynamics.NetworkDynamics.Graphs: path_graph
using PowerDynamics: Network, NWState, uflat, pflat
using PowerDynamics.NetworkDynamics: VIndex, EIndex, SII
using PowerDynamics.SciMLBase: SciMLBase, solve, ODEProblem, auto_dt_reset!
using PowerDynamics.NetworkDynamics.DiffEqCallbacks: PresetTimeCallback
using PowerDynamics.ModelingToolkit: @named

using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
using Makie: Makie, Figure, Axis, axislegend, lines!, Cycled

using PowerDynamics: PowerDynamics, compile_bus, compile_line, MTKLine
using PowerDynamics.Library: PiLine, SlackDifferential

using JLD2: JLD2
using Test: Test, @test, @test_broken

export TrajectoriesOfInterest, plottoi, compare
include("TrajectoriesOfInterest.jl")

export line_between_slacks, bus_on_slack
include("scenarios.jl")

export @reftest, set_reference_dir, refup
include("reftests.jl")

using PowerDynamics.Library: SauerPaiMachine, ConstantYLoad, AVRTypeI, TGOV1
using PowerDynamics: @initformula, CompositeInjector, MTKBus, pfSlack, pfPV, pfPQ, add_initformula!
include("testsystems.jl")


end # module PowerDynamicsTesting
