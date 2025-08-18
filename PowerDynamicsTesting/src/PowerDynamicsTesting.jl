module PowerDynamicsTesting

using Graphs: path_graph
using NetworkDynamics: Network, NWState, uflat, pflat
using NetworkDynamics: VIndex, EIndex, SII
using SciMLBase: SciMLBase, ODEProblem, solve, auto_dt_reset!
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks: PresetTimeCallback
using ModelingToolkit: @named
using Makie: Makie, Figure, Axis, axislegend, lines!, Cycled
using OrderedCollections: OrderedDict

using PowerDynamics: PowerDynamics, Bus, Line, MTKLine
using PowerDynamics.Library: PiLine, SlackDifferential

using JLD2: JLD2
using Test: Test, @test, @test_broken

export TrajectoriesOfInterest, plottoi, compare
include("TrajectoriesOfInterest.jl")

export line_between_slacks, bus_on_slack
include("scenarios.jl")

export @reftest, set_reference_dir, refup
include("reftests.jl")


end # module PowerDynamicsTesting
