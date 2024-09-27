module OpPoDyn

using NetworkDynamics: VertexFunction, EdgeFunction, ODEVertex, StaticEdge, Fiducial
using NetworkDynamics: has_default, set_default!, get_default
using ModelingToolkit: ModelingToolkit, ODESystem, structural_simplify, get_name
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

end
