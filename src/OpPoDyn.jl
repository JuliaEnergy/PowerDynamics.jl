module OpPoDyn

using NetworkDynamics: VertexFunction, EdgeFunction, ODEVertex, StaticEdge, Fiducial
using ModelingToolkit: ModelingToolkit, ODESystem, structural_simplify
using Graphs: SimpleGraph, add_edge!, nv
using ArgCheck: @argcheck
using Setfield: @set, @set!

export @attach_metadata!
include("utils.jl")
include("Library/Library.jl")
using .Library

export Line, Bus
export simplify_mtkline, simplify_mtkbus
include("network_components.jl")

end
