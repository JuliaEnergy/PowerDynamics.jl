module OpPoDynTesting

using Graphs: SimpleGraph, add_edge!, path_graph
using NetworkDynamics: Network, NWState, uflat, pflat
using NetworkDynamics: SymbolicIndex, VIndex, EIndex, vidxs, eidxs
using SciMLBase: SciMLBase, ODEProblem, solve, auto_dt_reset!
using OrdinaryDiffEqRosenbrock: Rodas5P, Rosenbrock23
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqNonlinearSolve: OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks: PresetTimeCallback
using ModelingToolkit: @named
using Makie: Makie, lines, Figure, Axis, axislegend, lines!, Cycled
using OrderedCollections: OrderedDict

using OpPoDyn: OpPoDyn, Bus, Line
using OpPoDyn.Library: DynawoPiLine, MTKLine, SlackDifferential

using JLD2: JLD2
using Test: Test, @test, @test_broken

export TrajectoriesOfInterest, plottoi, compare
include("TrajectoriesOfInterest.jl")

export line_between_slacks, bus_on_slack
include("scenarios.jl")

export @reftest, set_reference_dir, refup
include("reftests.jl")


end # module OpPoDynTesting
