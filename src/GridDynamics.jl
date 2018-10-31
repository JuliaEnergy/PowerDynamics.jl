# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using Base.Iterators: flatten

"Abstract super type for all abstract grid dynamics types."
abstract type GridDynamics end

"Abstract super type for all grid dynamics represented by ODEs."
abstract type AbstractOrdinaryGridDynamics <: GridDynamics end

"Abstract super type for all grid dynamics represented by DAEs."
abstract type AbstractAlgebraicGridDynamics <: GridDynamics end

NetworkRHS(g::GridDynamics) = g.rhs

"""
    struct OrdinaryGridDynamics <: AbstractOrdinaryGridDynamics
        rhs::NetworkRHS
    end

The data structure that contains all the information necessary for a power grid model that can be described
as an ordinary differential equation. In this case, only the [`PowerDynBase.NetworkRHS`] is necessary.
"""
@with_kw struct OrdinaryGridDynamics <: AbstractOrdinaryGridDynamics
    rhs::NetworkRHS
end
(dyn::OrdinaryGridDynamics)(dx_dt::AbstractVector, x_in::AbstractVector, p, t) = dyn.rhs(ODEVariable(x_in, dx_dt), t)

"""
    struct OrdinaryGridDynamicsWithMass <: AbstractAlgebraicGridDynamics
        rhs::NetworkRHS
        masses::AbstractVector{Bool} # diagonal part of the mass matrix, off-diagonal is assumed to be 0 anyway
    end

The data structure that contains all the information necessary for a power grid model that can be described
as an ordinary differential equation with masses, i.e. a semi-explicit differential algebraic equation.
`rhs` is the [`PowerDynBase.NetworkRHS`]. `masses` is a 1-dimensional array representing the diagonal entries
of the mass matrix. The off-diagonal entries are assumed to be 0. `masses` can only contain boolean values
representing: `true` the equation is treated as a ordinary differential eqation and `false` the equation is
treated as an algebraic constraint on the state variables.
"""
@with_kw struct OrdinaryGridDynamicsWithMass <: AbstractAlgebraicGridDynamics
    rhs::NetworkRHS
    masses::AbstractVector{Bool} # diagonal part of the mass matrix, off-diagonal is assumed to be 0 anyway
end
(dyn::OrdinaryGridDynamicsWithMass)(dx_dt::AbstractVector, x_in::AbstractVector, p, t) = dyn.rhs(ODEVariable(x_in, dx_dt), t)

"""
    struct AlgebraicGridDynamics <: AbstractAlgebraicGridDynamics
        rhs::NetworkRHS
        differentials::AbstractVector{Bool} # boolean values whether there a variable is a differential
    end

The data structure that contains all the information necessary for a power grid model that can be described
as an differential algebraic equation.
`rhs` is the [`PowerDynBase.NetworkRHS`]. `differentials` is a 1-dimensional array of boolean values.
A `true` entry means the corresponding variable is dynamic and has a derivative variable.
A `false` entry means the corresponding variable is defined by an algebraic constraint only.
"""
@with_kw struct AlgebraicGridDynamics <: AbstractAlgebraicGridDynamics
    rhs::NetworkRHS
    differentials::AbstractVector{Bool} # boolean values whether there a variable is a differential
end
(dyn::AlgebraicGridDynamics)(x_out::AbstractVector, dx_dt::AbstractVector, x_in::AbstractVector, p, t) = dyn.root(DAEVariable(val=x_in, ddt=dx_dt, out=x_out), t)


"Create for each subtype of [`PowerDynBase.AbstractNodeDynamics`](#ref) the corresponding subtype of [`PowerDynBase.GridDynamics`](#ref)."
_GridDynamics(nodes::AbstractVector{<:OrdinaryNodeDynamics}, LY::AbstractMatrix) =
    OrdinaryGridDynamics(NetworkRHS(nodes, LY))

_GridDynamics(nodes::AbstractVector{<:OrdinaryNodeDynamicsWithMass}, LY::AbstractMatrix) =
    OrdinaryGridDynamicsWithMass(NetworkRHS(nodes, LY), masses(nodes))

_GridDynamics(node_dynamics::AbstractVector{<:AbstractNodeDynamics}, ::AbstractMatrix) =
    throw(GridDynamicsError("There is seemingly no grid dynamics constructor implemented for your combination of node dynamics (or there is a bug)."))

"""Check whether the admittance laplacian has no purely nodal admittances,
i.e. that the sum of columns and rows equals to zero."""
function checkLY(LY::AbstractMatrix)
    if ~( sum(LY, dims=[1]) .|> approxzero |> all)
        throw(GridDynamicsError("rows of the admittance matrix have to sum to 0"))
    end
    if ~( sum(LY, dims=[2]) .|> approxzero |> all)
        throw(GridDynamicsError("columns of the admittance matrix have to sum to 0"))
    end
end

"""
    function GridDynamics(
        nodes::AbstractArray{<:AbstractNodeDynamics},
        LY::AbstractMatrix;
        skip_LY_check=false,
        kwargs...)

Bring all sutypes of [`PowerDynBase.AbstractNodeDynamics`](#ref) on one level and then
create for each subtype of [`PowerDynBase.AbstractNodeDynamics`](#ref) the corresponding subtype of [`PowerDynBase.GridDynamics`](#ref) by using
[`PowerDynBase._GridDynamics`](#ref).
"""
function GridDynamics(
    nodes::AbstractArray{<:AbstractNodeDynamics},
    LY::AbstractMatrix;
    skip_LY_check=false,
    kwargs...
    )
    if !skip_LY_check
        checkLY(LY)
    end
    common_type = promote_type(map(typeof, nodes)...) # get the common node dynamics type
    _GridDynamics(convert(Array{common_type}, nodes), LY; kwargs...)
end


"""
    function GridDynamics(nodes::AbstractVector{<:AbstractNodeParameters}, args...; kwargs...)

Convert all subtypes of [`PowerDynBase.AbstractNodeParameters`](#ref) to the corresponding
subtypes of [`PowerDynBase.AbstractNodeDynamics`](#ref) and then call [`PowerDynBase.GridDynamics`](#ref)
again.
"""
function GridDynamics(nodes::AbstractVector{<:AbstractNodeParameters}, args...; kwargs...)
    GridDynamics(map(construct_node_dynamics, nodes), args...; kwargs...)
end

"Return the 1-dimensional differentials array (see [`PowerDynBase.AlgebraicGridDynamics`]) for the internal variables for each node."
int_differentials(node::AbstractAlgebraicNodeDynamics, nodes::AbstractAlgebraicNodeDynamics...) = [int_differentials(node); int_differentials(nodes...)]
"Return the 1-dimensional differentials array (see [`PowerDynBase.AlgebraicGridDynamics`]) for the voltagee variables for each node."
u_differentials(node::AbstractAlgebraicNodeDynamics, nodes::AbstractAlgebraicNodeDynamics...) = [u_differentials(node); u_differentials(nodes...)]
u_differentials(node::OrdinaryNodeDynamicsWithMass) = repeat([node.m_u], inner=2) # 1 Complex number is 2 Real numbers
int_differentials(node::OrdinaryNodeDynamicsWithMass) = node.m_int
u_differentials(node::AlgebraicNodeDynamics) = repeat([node.d_u], inner=2)
int_differentials(node::AlgebraicNodeDynamics) = node.d_int

"Collect the differnetials (or masses) for the whole grid dynamics from each node dynamics."
function differentials(nodes::AbstractVector{<:AbstractAlgebraicNodeDynamics})
    return [u_differentials(nodes...) ; int_differentials(nodes...)]
end

differentials(rhs::NetworkRHS) = differentials(Nodes(rhs))
differentials(g::OrdinaryGridDynamicsWithMass) = differentials(NetworkRHS(g))

const masses = differentials # structurally the same, only the interpretation is different
