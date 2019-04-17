# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: convert, promote_rule

"Abstract super type for all abstract node dynamics types."
abstract type AbstractNodeDynamics{N <: AbstractNodeParameters} end

"Abstract super type for all node dynamics represented by ODEs."
abstract type AbstractOrdinaryNodeDynamics{N <: AbstractNodeParameters} <: AbstractNodeDynamics{N} end

"Abstract super type for all node dynamics represented by DAEs."
abstract type AbstractAlgebraicNodeDynamics{N <: AbstractNodeParameters} <: AbstractNodeDynamics{N} end

function MassMatrix(m_u::Bool = false, m_int = no_internal_masses)
    mm = [m_u, m_u] # double mass for u, because it is a complex variable
    append!(mm, m_int)
    return mm .|> Int
end

################################################################################
# ODEs
################################################################################

@doc doc"""
```Julia
OrdinaryNodeDynamics(;rhs, n_int)
```

The type representing the dynamics of a node that is described via ODEs.

Each node ``a`` has the complex voltage ``u`` and ``n`` real internal variables ``y_1, \dots, y_n``, so it
generally describes a system of ordinary differential equation as

```math
\frac{du_a}{dt} = f_u(u_a, {i_c}_a, y_1, \dots, y_n) \\
\frac{dy_{ak}}{dt} = f_k(u_a, {i_c}_a, y_1, \dots, y_n)\quad \forall k = 1, \dots, n.
```
``f`` is represented by `rhs` field of `OrdinaryNodeDynamics`.
- the general signature of `rhs` is
```Julia
rhs(dint_dt::AbstractVector,
    u::Complex,
    i::Complex,
    int::AbstractVector,
    t,
    )::Complex
```
- Input
  - `u` is the complex voltage ``u``
  - `i` is the complex current ``i``
  - `int` is the array of internal variables ``y_1, \dots, y_n``
  - `t` is the time ``t``
- Output
  - the (complex) return value describes ``\frac{du}{dt}``
  - `rhs` writes values in `dint_dt` describing the left-hand side ``\frac{dy_1}{dt}, \dots, \frac{dy_n}{dt}``

"""
@with_kw struct OrdinaryNodeDynamics{N <: AbstractNodeParameters} <: AbstractOrdinaryNodeDynamics{N}
    rhs::Function # how to define the function type, should be clear so the interface is forced, keyword FunctionWrapper
    parameters::N
    n_int
end
function (dyn::OrdinaryNodeDynamics)(n, u::ODEVariable, i,
    int::ODEVariable , t)
    u.ddt[n] =  dyn.rhs(int.ddt, u.val[n], i[n], int.val, t)
    nothing
end

"Identify each subtype of [`AbstractNodeDynamics`](@ref) with its corresponding subtype of [`AbstractDEVariable`](@ref)"
getDEVariableType(::Type{Val{OrdinaryNodeDynamics}}) = ODEVariable

"Get number of internal arguments of the node."
nint(dyn::OrdinaryNodeDynamics) = dyn.n_int

"Get the parameters struct for the node."
parametersof(n::OrdinaryNodeDynamics) = n.parameters

################################################################################
# ODEs with mass matrix
################################################################################

"A variable to be used when no internal masses are present for a node dynamics type."
const no_internal_masses = Vector{Bool}()

@doc doc"""
```Julia
OrdinaryNodeDynamicsWithMass(;rhs, n_int, m_u, m_int)
```

The type representing the dynamics of a node that is described via ODEs.

Each node ``a`` has the complex voltage ``u`` and ``n`` (`= n_int`) real internal variables ``y_1, \dots, y_n``, so it
generally describes a system of ordinary differential equation with
a voltage mass ``m_u`` and internal masses ``m^{int}_1, \dots, m^{int}_n`` as

```math
m_u\frac{du_a}{dt} = f_u(u_a, {i_c}_a, y_1, \dots, y_n) \\
m^{int}_k\frac{dy_{ak}}{dt} = f_k(u_a, {i_c}_a, y_1, \dots, y_n)\quad \forall k = 1, \dots, n.
```

As we assume that all masses are binary (either 1, or 0), that means, one can implement [semi-explicit differential algebraic equations](https://en.wikipedia.org/wiki/Differential-algebraic_system_of_equations) with
this node dynamics type.
``f`` is represented by `rhs` field of `OrdinaryNodeDynamics`.
- the general signature of `rhs` is
```Julia
rhs(dint_dt::AbstractVector,
    u::Complex,
    i::Complex,
    int::AbstractVector,
    t,
    )::Complex
```
- Input
  - `u` is the complex voltage ``u``
  - `i` is the complex current ``i``
  - `int` is the array of internal variables ``y_1, \dots, y_n``
  - `t` is the time ``t``
- Output
  - the (complex) return value describes ``\frac{du}{dt}``
  - `rhs` writes values in `dint_dt` describing the left-hand side ``\frac{dy_1}{dt}, \dots, \frac{dy_n}{dt}``

The binary masses are:
- `m_u` is the boolean value for ``m_u``
- `m_int` is the array of boolean values for ``m^{int}_1, \dots, m^{int}_n``

"""
struct OrdinaryNodeDynamicsWithMass{N <: AbstractNodeParameters} <: AbstractAlgebraicNodeDynamics{N}
    ode_dynamics::OrdinaryNodeDynamics{N}
    m_u::Bool # Answers the question: Is the voltage treated as a dynamic variable with a differential
    m_int::AbstractVector{Bool} # for each internal variable: true if there is a differential for it, else false (if it is an algebraic constraint only)
end
OrdinaryNodeDynamicsWithMass(;rhs::Function, n_int, m_u, m_int, parameters) = OrdinaryNodeDynamicsWithMass(
    OrdinaryNodeDynamics(rhs, parameters, n_int),
    m_u, m_int)
(dyn::OrdinaryNodeDynamicsWithMass)(args...;kwargs...) = dyn.ode_dynamics(args...;kwargs...)

getDEVariableType(::Type{Val{OrdinaryNodeDynamicsWithMass}}) = ODEVariable

nint(dyn::OrdinaryNodeDynamicsWithMass) = nint(dyn.ode_dynamics)

OrdinaryNodeDynamics(n::OrdinaryNodeDynamicsWithMass) = n.ode_dynamics

# TODO: remove OrdinaryNodeDynamics(::OrdinaryNodeDynamicsWithMass) and access the field here directly
parametersof(n::OrdinaryNodeDynamicsWithMass) = n |> OrdinaryNodeDynamics |> parametersof

"""
    convert(::Type{OrdinaryNodeDynamicsWithMass}, ::OrdinaryNodeDynamics)


Conversion of `OrdinaryNodeDynamics` to `OrdinaryNodeDynamicsWithMass` by assuming all masses are 1 (`true`).
"""
convert(::Type{OrdinaryNodeDynamicsWithMass}, dyn::OrdinaryNodeDynamics) =
    OrdinaryNodeDynamicsWithMass(dyn, true, ones(Bool, dyn.n_int))

"""
    promote_rule(::Type{OrdinaryNodeDynamics}, ::Type{OrdinaryNodeDynamicsWithMass}) = OrdinaryNodeDynamicsWithMass

`OrdinaryNodeDynamics` can be promoted to `OrdinaryNodeDynamicsWithMass`, see [`PowerDynBase.convert`](@ref).
"""
promote_rule(::Type{<:OrdinaryNodeDynamics}, ::Type{<:OrdinaryNodeDynamicsWithMass}) =
    OrdinaryNodeDynamicsWithMass


################################################################################
# DAEs
################################################################################

"A variable to be used when no internal differentials are present for a node dynamics type."
const no_internal_differentials = Vector{Bool}()

"DOCS TBD!"
@with_kw struct AlgebraicNodeDynamics{N <: AbstractNodeParameters} <: AbstractAlgebraicNodeDynamics{N}
    root::Function # how to define the function type, should be clear so the interface is forced, keyword FunctionWrapper
    parameters::N
    n_int
    d_u::Bool # Answers the question: Is the voltage treated as a dynamic variable with a differential
    d_int::AbstractVector{Bool} # for each internal variable: true if there is a differential for it, else false (if it is an algebraic constraint only)
end
function (dyn::AlgebraicNodeDynamics)(i, u::DAEVariable, i_c,
    int::DAEVariable , t)
    u.out[i] =  dyn.root(int.out, u.ddt[i], int.ddt, u.val[i], i_c[i], int.val, t)
    nothing
end

getDEVariableType(::Type{Val{AlgebraicNodeDynamics}}) = DAEVariable

nint(dyn::AlgebraicNodeDynamics) = dyn.n_int

parametersof(n::AlgebraicNodeDynamics) = n.parameters

"A function converting a rhs-type function to a root-type function."
function rhs2root(rhs::Function)
    function root_from_rhs!(
        int_out::AbstractVector,
        du, #::Complex
        dint::AbstractVector,
        u, #::Complex
        i_c, #::Complex
        int::AbstractVector,
        t,
        )#::Complex
        rhs_u = rhs(int_out, u, i_c, int, t)
        int_out[:] -= dint
        return rhs_u - du
      end
end

"""
    convert(::Type{AlgebraicNodeDynamics}, ::OrdinaryNodeDynamics)

Conversion of `OrdinaryNodeDynamics` to `AlgebraicNodeDynamics` by going via `OrdinaryNodeDynamicsWithMass`.
"""
convert(::Type{AlgebraicNodeDynamics}, node_dynamics::OrdinaryNodeDynamics) = @>> node_dynamics convert(OrdinaryNodeDynamicsWithMass) convert(AlgebraicNodeDynamics)

"""
    convert(::Type{AlgebraicNodeDynamics}, ::OrdinaryNodeDynamicsWithMass)

Conversion of `OrdinaryNodeDynamicsWithMass` to `AlgebraicNodeDynamics` by transforming a RHS function into a root function.
"""
function convert(::Type{AlgebraicNodeDynamics}, node::OrdinaryNodeDynamicsWithMass)
    AlgebraicNodeDynamics(
        root=rhs2root(node.ode_dynamics.rhs),
        parameters=parametersof(node),
        n_int=node.ode_dynamics.n_int,
        d_u=node.m_u,
        d_int=node.m_int
    )
end

"""
    promote_rule(::Type{OrdinaryNodeDynamics}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics

`OrdinaryNodeDynamics` can be promoted to `AlgebraicNodeDynamics`, see [`PowerDynBase.convert`](@ref).
"""
promote_rule(::Type{OrdinaryNodeDynamics}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics

"""
    promote_rule(::Type{OrdinaryNodeDynamicsWithMass}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics

`OrdinaryNodeDynamicsWithMass` can be promoted to `AlgebraicNodeDynamics`, see [`PowerDynBase.convert`](@ref).
"""
promote_rule(::Type{OrdinaryNodeDynamicsWithMass}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics
