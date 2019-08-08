# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: @__doc__
using Base.Meta
using Base.Iterators
using MacroTools
using NetworkDynamics
using LinearAlgebra

function cndfunction_builder!(
    internals,
    massmatrix,
    func_body,
    cndfunction
    )
    rhscall = :(rhs!(
        dx,
        x,
        e_s,
        e_d,
        p,
        t)
        )
    rhsbody = quote end
    rhsbody.args[1] = func_body.args[1]
    append!(rhsbody.args, [:(i = total_current(e_s, e_d))])
    append!(rhsbody.args, [:(u = x[1] + x[2] * im)])

    extracted_vars =  [:($sym = x[$(index+2)] ) for (index, sym) in enumerate(internals.vars)]
    append!(rhsbody.args, extracted_vars)
    append!(rhsbody.args, func_body.args)
    extracted_dvars =  [:(dx[$(index+2)] = $sym) for (index, sym) in enumerate(internals.dvars)]
    du_real = [:(dx[1] = real(du))]
    du_imag = [:(dx[2] = imag(du))]
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
    append!(catchdef.args[1].args, [du_real; du_imag ; extracted_dvars; :(return nothing)])
    append!(rhsbody.args, [catchdef])
    rhsfunction = Expr(:function, rhscall, rhsbody)
    all_syms = [:u_r, :u_i]
    append!(all_syms, internals.vars)
    ode_vertex = :(ODEVertex(f! = rhs!, dim = $(length(internals.vars) +2), mass_matrix=$massmatrix, sym = $all_syms))
    append!(cndfunction.args[2].args, [rhsfunction, ode_vertex])

    nothing
end


function buildparameterstruct(name, parameters)
    struct_def = Expr(
        :struct, false,
        :($name), # define the struct as a subtype of AbstractNodeParameters
        Expr(:block, parameters..., # set all the parmeters as fields in the struct
            Expr(:(=), # define the constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...)
            )
        )
    )
end

function getinternalvars(internalsdef)
    (vars = map(ex -> ex.args[1] , internalsdef.args),
        dvars = map(ex -> ex.args[2] , internalsdef.args))
end

function generate_symbolsof_fct(name, internals)
    vars = [:u_r, :u_i]
    append!(vars, internals.vars)
    :(symbolsof(::$name) = $(Expr(:vect, QuoteNode.(vars)...)))
end

function generate_dimension_fct(name, internals)
    :(dimension(::$name) = $(length(internals.vars)+2))
end

"""See [`PowerDynamics.@DynamicNode`](@ref)."""
function DynamicNode(typedef, massmatrix, prep, internalsdef, func_body)
    @capture(typedef, name_(parameters__))
    internals = getinternalvars(internalsdef)
    # build parameters struct
    struct_def = buildparameterstruct(name, parameters)

    massmatrix = massmatrix === nothing ? I : massmatrix

    # build `construct_vertex`
    cndcall = :(construct_vertex(par::$(name)))
    extracted_parameters = map(sym -> :( $sym = par.$sym ), parameters)
    cndbody = quote end
    append!(cndbody.args, extracted_parameters)
    append!(cndbody.args, prep.args)
    cndfunction = Expr(:function, cndcall, cndbody)

    # this might be different for
    cndfunction_builder!(
        internals,
        massmatrix,
        func_body,
        cndfunction
        )

    fct_symbolsof = generate_symbolsof_fct(name, internals)
    fct_dimension = generate_dimension_fct(name, internals)

    ret = quote
        @__doc__ $(struct_def)
        $(cndfunction)
        $(fct_symbolsof)
        $(fct_dimension)
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
    return create(typedef, nothing, prep, internals, func_body)
end

macro DynamicNode(typedef, massmatrix_exp, prep, internals, func_body)
    massmatrix = eval(massmatrix_exp)
    print(massmatrix_exp)
    return create(typedef, massmatrix, prep, internals, func_body)
end

function create(typedef, massmatrix, prep, internals, func_body)
    mainex = DynamicNode(typedef, massmatrix, prep, internals, func_body)
    @capture(typedef, name_(parameters__))
    mainexstr = "$(copy(mainex)|>rmlines|> MacroTools.striplines)"
    showex = :(showdefinition(io::IO, ::Type{$name}) = println(io, $mainexstr))
    append!(mainex.args, [showex])
    return esc(mainex)
end

"Show the definition generated by a macro of PowerDynamics.jl,
e.g. the macro [`PowerDynamics.@DynamicNode`](@ref)"
showdefinition(x) = showdefinition(stdout, x)

"A variable to be used when no internal masses are present for a node dynamics type."
const no_internal_masses = Vector{Bool}()

function MassMatrix(;m_u::Bool = false, m_int = no_internal_masses)
    mm = [m_u, m_u] # double mass for u, because it is a complex variable
    append!(mm, m_int)
    return mm .|> Int
end
