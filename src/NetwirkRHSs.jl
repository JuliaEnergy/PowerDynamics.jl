# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

"Get the total number of internal variables for an array of node dynamics."
total_nint(nodes::AbstractVector{<:AbstractNodeDynamics}) = sum(nint, nodes)

"""
Get the total number of dynamic variables for an array of node dynamics.

This is basically the (real) dimension of the system, hence the sum of internal
dynamic variables + 2*(number of nodes = number of complex voltages). The 2 is
due to the fact that the complex voltages are treated as two real variables.
"""
total_nvars(nodes::AbstractVector{<:AbstractNodeDynamics}) = sum(nint, nodes) + 2*length(nodes)

"Get the unit ranges that indicate where in the array the internal variables for each of the nodes is saved."
function internal_unitranges(nodes::AbstractVector{<:AbstractNodeDynamics})
    next_unit_range(range, node) = (range.stop+1):(range.stop+nint(node))
    return accumulate(next_unit_range, nodes, init=1:0)
end

"""
    abstract type AbstractNetworkFunction{T<:AbstractNodeDynamics, M<:AbstractMatrix} end

Abstract super type of all functions that define how a differential equation for the whole
network / power grid behaves, e.g. the full right-hand-side function of the ODE.
"""
abstract type AbstractNetworkFunction{T<:AbstractNodeDynamics, M<:AbstractMatrix} end
"""
    struct NetworkRHS{T, M} <: AbstractNetworkFunction{T, M}
        nodes::AbstractVector{T}
        LY::M
        numnodes
        systemsize
        intrange # unitrange telling me where I find the internal dynamic variables
        nodalintranges # unit ranges to find the internal variables for each node in the full length of internal variables
        rhsinterface::Function
    end

Representing the full dynamics of the power grid.

From [DPSABase.`AbstractNetworkFunction`](#ref) it inherits the conditions that `T<:AbstractNodeDynamics` and `M<:AbstractMatrix}`.
"""
@with_kw struct NetworkRHS{T, M} <: AbstractNetworkFunction{T, M}
    nodes::AbstractVector{T}
    LY::M
    numnodes
    systemsize
    intrange # unitrange telling me where I find the internal dynamic variables
    nodalintranges # unit ranges to find the internal variables for each node in the full length of internal variables
end
function (rhs::NetworkRHS)(x::AbstractDEVariable, t)
    # distribute the values of the of the DEVariables over all
    nodeiterator(rhs, x, t)
    nothing
end
# external constructor
function NetworkRHS(nodes::AbstractVector{T}, LY::M) where {T<:AbstractNodeDynamics, M<:AbstractMatrix}
    numnodes = length(nodes)
    @assert size(LY) == (numnodes, numnodes)
    systemsize = total_nvars(nodes)
    NetworkRHS(
        nodes = nodes,
        LY = LY,
        numnodes = numnodes,
        systemsize = systemsize,
        intrange = (2 * numnodes + 1):systemsize,
        nodalintranges = internal_unitranges(nodes)
    )
end
Nodes(rhs::NetworkRHS) = rhs.nodes
SystemSize(rhs::NetworkRHS) = rhs.systemsize
AdmittanceLaplacian(rhs::NetworkRHS) = rhs.LY

# fall-back: only GridDynamics(x) needs to be defined, the rest is automatically via this line and the following
NetworkRHS(x) = x |> GridDynamics |> NetworkRHS
SystemSize(x) = x |> NetworkRHS |> SystemSize
Nodes(x) = x |> NetworkRHS |> Nodes
AdmittanceLaplacian(x) = x |> NetworkRHS |> AdmittanceLaplacian

"""
    nodeiterator(rhs::NetworkRHS, x::AbstractDEVariable, t)

Distribute the values in `x` over all the nodes that are summarized in `rhs`.
"""
function nodeiterator(rhs::NetworkRHS, x::AbstractDEVariable, t)
    # node the double parenthesis as it gets a tuple because that makes the
    # concatenation easier
    @assert length(x) == rhs.systemsize " not( length(x)=$(length(x)) == system_size=$system_size)"
    u = complexview(x, 1, length(Nodes(rhs))) # complex voltages
    int = view(x, rhs.intrange)
    i = rhs.LY * u.val # complex current
    for (n, n_int_range) in enumerate(rhs.nodalintranges)
        Nodes(rhs)[n](n, u, i, view(int, n_int_range) , t)
    end
end
