module OpPoDyn

using NetworkDynamics: NetworkDynamics, VertexModel, EdgeModel,
                       has_default, set_default!, get_default,
                       has_graphelement, set_graphelement!, get_graphelement, Network,
                       has_metadata, get_metadata, set_metadata!, NWState, uflat, pflat
using ModelingToolkit: ModelingToolkit, ODESystem, structural_simplify, get_name, getname, @named
using SciMLBase: SciMLBase, solve
using NonlinearSolve: NonlinearSolve, NonlinearProblem
using ForwardDiff: ForwardDiff
using LinearAlgebra: LinearAlgebra, Diagonal, diag, pinv, eigvals
using Graphs: SimpleGraph, add_edge!, nv
using ArgCheck: @argcheck
using Setfield: @set, @set!

export @attach_metadata!, set_voltage!, set_current!
include("utils.jl")
include("Library/Library.jl")
using .Library

export Line, Bus
export simplify_mtkline, simplify_mtkbus
include("network_components.jl")

using DataFrames: DataFrame
using OrderedCollections: OrderedDict
export pfSlack, pfPV, pfPQ
include("powerflow.jl")

end
