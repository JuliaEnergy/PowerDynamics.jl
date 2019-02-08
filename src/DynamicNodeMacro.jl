# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: @__doc__
using Base.Meta
using Base.Iterators
using MacroTools

function cndfunction_builder!(::Type{Val{:OrdinaryNodeDynamics}},
    internals,
    dynamicscall,
    func_body,
    cndfunction
    )
    rhscall = :(rhs!(
        dint::AbstractVector, # output for internal variables
        u,#::Complex,
        i,#::Complex,
        int::AbstractVector,
        t)#::Complex
        )
    rhsbody = quote end
    rhsbody.args[1] = func_body.args[1]
    extracted_vars =  [:($sym = int[$index] ) for (index, sym) in enumerate(internals.vars)]
    append!(rhsbody.args, extracted_vars)
    append!(rhsbody.args, func_body.args)
    extracted_dvars =  [:(dint[$index] = $sym) for (index, sym) in enumerate(internals.dvars)]
    errorif = Expr(:block,
        :(if typeof(e) === UndefVarError
            throw(NodeDynamicsError("you need to provide $(e.var)"))
        else
            throw(e)
        end))
    # delete the line statement so the NodeDynamicsError points to the
    # actual node dynamics definition and not here (:
    # That way a user is not confused (hopefully)
    deleteat!(errorif.args[1].args[2].args, 1)
    catchdef = Expr(:try, Expr(:block), :e, errorif)
    append!(catchdef.args[1].args, [extracted_dvars; :(return du)])
    append!(rhsbody.args, [catchdef])
    rhsfunction = Expr(:function, rhscall, rhsbody)
    append!(dynamicscall.args, [
        Expr(:kw, :rhs, :rhs!),
        Expr(:kw, :n_int, length(internals.vars)),
        Expr(:kw, :parameters, :par),
        ])
    append!(cndfunction.args[2].args, [rhsfunction, dynamicscall])

    nothing
end

function cndfunction_builder!(::Type{Val{:OrdinaryNodeDynamicsWithMass}}, args...; kwargs...)
    cndfunction_builder!(Val{:OrdinaryNodeDynamics}, args...; kwargs...)
end

function cndfunction_builder!(::Type{Val{T}}, args...;kwargs...) where {T}
    throw(NodeDynamicsError("unknown node dynamics type $T"))
end

function buildparameterstruct(name, parameters)
    struct_def = Expr(
        :struct, false,
        :($name <: AbstractNodeParameters), # define the struct as a subtype of AbstractNodeParameters
        Expr(:block, parameters..., # set all the parmeters as fields in the struct
            Expr(:(=), # define the constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...)
            )
        )
    )
end

function getinternalvars(::Type{Val{:OrdinaryNodeDynamics}}, internalsdef)
    (vars = map(ex -> ex.args[1] , internalsdef.args),
        dvars = map(ex -> ex.args[2] , internalsdef.args))
end

function getinternalvars(::Type{Val{:OrdinaryNodeDynamicsWithMass}}, internalsdef)
    getinternalvars(Val{:OrdinaryNodeDynamics}, internalsdef)
end

function getinternalvars(::Type{Val{T}}, args...;kwargs...) where {T}
    throw(NodeDynamicsError("cannot extract internal symbols for $T. Are you sure that is implemented?"))
end

function generate_symbolsof_fct(::Type{Val{:OrdinaryNodeDynamics}}, name, internals)
    :(symbolsof(::Type{$name}) = ODENodeSymbols($(Expr(:vect, QuoteNode.(internals.vars)...)), $(Expr(:vect, QuoteNode.(internals.dvars)...))))
end

function generate_symbolsof_fct(::Type{Val{:OrdinaryNodeDynamicsWithMass}}, name, internals)
    generate_symbolsof_fct(Val{:OrdinaryNodeDynamics}, name, internals)
end

function generate_symbolsof_fct(::Type{Val{T}}, name, internals) where T
    throw(NodeDynamicsError("cannot generator `symbolsof` function for $T. Are you sure that is implemented?"))
end

"""See [`PowerDynBase.@DynamicNode`](@ref)."""
function DynamicNode(typedef, prep, internalsdef, func_body)
    @capture(typedef, name_(parameters__) <: dynamicscall_)
    dynamicstype = dynamicscall.args[1]
    internals = getinternalvars(Val{dynamicstype}, internalsdef)

    # build parameters struct
    struct_def = buildparameterstruct(name, parameters)

    # build `construct_node_dynamics`
    cndcall = :(construct_node_dynamics(par::$(name)))
    extracted_parameters = map(sym -> :( $sym = par.$sym ), parameters)
    cndbody = quote end
    append!(cndbody.args, extracted_parameters)
    append!(cndbody.args, prep.args)
    cndfunction = Expr(:function, cndcall, cndbody)

    # this might be different for
    cndfunction_builder!(Val{dynamicstype},
        internals,
        dynamicscall,
        func_body,
        cndfunction
        )

    fct_symbolsof = generate_symbolsof_fct(Val{dynamicstype}, name, internals)

    ret = quote
        @__doc__ $(struct_def)
        $(cndfunction)
        $(fct_symbolsof)
    end
    return ret
end

"""
Macro for creating a new type of dynmic nodes.

Syntax Description:
```Julia
@DynamicNode MyNewNodeName(Par1, Par2, ...) <: NodeDynamicsType(N1, N2, ...) begin
    [all prepratory things that need to be run just once]
end [[x1, dx1], [x2, dx2]] begin
    [the actual dynamics equation]
    [important to set the output variables]
end
```
where `MyNewNodeName` is the name of the new type of dynamic node, `Par1, Par2, ...`
are the names of the parameters, `NodeDynamicsType` the the node dynamics type
(e.g. `OrdinaryNodeDynamics` or `OrdinaryNodeDynamicsWithMass`), `N1, N1, ...`
the parameters of the dynamics type, `x1, x2, ...` the internal variables of the
node and `dx1, dx2, ...` the corresponding differentials.

In the first block, the preparation code that needs to be run only once is inserted.
Finally, the second block contains the dynamics description, where it's important
that the output variables need to be set. In case of `OrdinaryNodeDynamics` and
`OrdinaryNodeDynamicsWithMass`, these are `du` and the differentials of the
internal variables (here `dx1, dx2`).

Below are two examples:

```Julia
@DynamicNode SwingEqParameters(H, P, D, Ω) <: OrdinaryNodeDynamics() begin
    @assert D > 0 "damping (D) should be >0"
    @assert H > 0 "inertia (H) should be >0"
    Ω_H = Ω * 2pi / H
end [[ω, dω]] begin
    p = real(u * conj(i_c))
    dϕ = ω # dϕ is only a temp variable that Julia should optimize out
    du = u * im * dϕ
    dω = (P - D*ω - p)*Ω_H
end
```

```Julia
@DynamicNode SlackAlgebraicParameters(U) <: OrdinaryNodeDynamicsWithMass(m_u=false, m_int=no_internal_masses) begin
end [] begin
        du = u - U
end
```

"""
macro DynamicNode(typedef, prep, internals, func_body)
    mainex = DynamicNode(typedef, prep, internals, func_body)
    @capture(typedef, name_(parameters__) <: dynamicscall_)
    mainexstr = "$(copy(mainex)|>rlr)"
    showex = :(showdefinition(io::IO, ::Type{$name}) = println(io, $mainexstr))
    append!(mainex.args, [showex])
    return esc(mainex)
end

"Show the definition generated by a macro of PowerDynamics.jl,
e.g. the macro [`PowerDynBase.@DynamicNode`](@ref) creating a subtype of [`PowerDynBase.AbstractNodeParameters`](@ref)."
showdefinition(x) = showdefinition(stdout, x)
