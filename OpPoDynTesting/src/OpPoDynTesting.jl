module OpPoDynTesting

using Graphs: SimpleGraph, add_edge!, path_graph
using NetworkDynamics: Network, NWState, uflat, pflat
using NetworkDynamics: SymbolicIndex, VIndex, EIndex, vidxs, eidxs
using SciMLBase: SciMLBase, ODEProblem, solve, auto_dt_reste!
using OrdinaryDiffEqRosenbrock: Rodas5P, Rosenbrock23
using OrdinaryDiffEqTsit5: Tsit5
using DiffEqCallbacks: PresetTimeCallback
using ModelingToolkit: @named
using Makie: lines, Figure, Axis, axislegend, lines!, Cycled
using OrderedCollections: OrderedDict

using OpPoDyn
using OpPoDyn.Library

export TrajectoriesOfInterest, plottoi, similartoi
include("utils.jl")

export line_between_slacks, bus_on_slack
include("scenarios.jl")

end # module OpPoDynTesting
