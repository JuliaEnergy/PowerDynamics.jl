# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using Base.Meta
using Base.Iterators
using MacroTools
using NetworkDynamics
using LinearAlgebra

function cndfunction_builder!(
    internals,
    massmatrix,
    func_body,
    cndfunction)
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
    append!(rhsbody.args, [:(i = total_current(e_s, e_d)+Y_n*(x[1] + x[2] * im))])
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
        :($name <: AbstractNode),
        Expr(:block, parameters...) # set all the parmeters as fields in the struct
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
    full_params = push!(copy(parameters), :Y_n)
    struct_def = buildparameterstruct(name, full_params)

    kw_constructor = Expr(:(=), # define the constructor
        Expr(:call, name, Expr(:parameters, parameters..., Expr(:kw, :Y_n, 0))),
        Expr(:call, name,  (full_params)...)
    )

    massmatrix = massmatrix === nothing ? I : massmatrix

    # build `construct_vertex`
    cndcall = :(construct_vertex(par::$(name)))
    extracted_parameters = map(sym -> :( $sym = par.$sym ), full_params)
    cndbody = quote end
    append!(cndbody.args, extracted_parameters)
    append!(cndbody.args, prep.args)
    cndfunction = Expr(:function, cndcall, cndbody)

    # this might be different for
    cndfunction_builder!(
        internals,
        massmatrix,
        func_body,
        cndfunction)

    fct_symbolsof = generate_symbolsof_fct(name, internals)
    fct_dimension = generate_dimension_fct(name, internals)

    ret = quote
        Base.@__doc__($(struct_def))
        $(kw_constructor)
        $(cndfunction)
        $(fct_symbolsof)
        $(fct_dimension)
    end
    return ret
end

"""
Macro for creating a new type of dynamic nodes.

Syntax Description:
```Julia
@DynamicNode MyNewNodeName(Par1, Par2, ...) begin
    [MassMatrix definition]
end begin
    [all prepratory things that need to be run just once]
end [[x1, dx1], [x2, dx2]] begin
    [the actual dynamics equation]
    [important to set the output variables]
end
```
where `MyNewNodeName` is the name of the new type of dynamic node, `Par1, Par2, ...`
are the names of the parameters, `x1, x2, ...` the internal variables of the
node and `dx1, dx2, ...` the corresponding differentials.

In the first block a MassMatrix may be specified. Using the [`MassMatrix`](@ref)
helper function here is recommended.
The whole block can be omitted and the identity matrix I is then used as default.

In the second block, the preparation code that needs to be run only once is inserted.
Finally, the third block contains the dynamics description, where it's important
that the output variables need to be set. These are `du` and the differentials of the
internal variables (here `dx1, dx2`).

Below are two examples:

```Julia
@DynamicNode SwingEqParameters(H, P, D, Ω) begin
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
@DynamicNode SlackAlgebraicParameters(U) begin
    MassMatrix() # no masses
end begin
    # empty prep block
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

"""
    MassMatrix(;m_u::Bool = false, m_int = no_internal_masses)

Creates a massmatrix.
Calling:
```Julia
MassMatrix()
```
creates a mass matrix with all masses turned off.

# Keyword Arguments
- `m_u::Bool=false`: Mass matrix value for the complex voltage u.
- `m_int::Vector{Bool}=Vector{Bool}()`: A vector representing the diagonal of the mass matrix. Specifies the masses for the remaining/internal variables.

"""
function MassMatrix(;m_u::Bool = false, m_int = no_internal_masses)
    mm = [m_u, m_u] # double mass for u, because it is a complex variable
    append!(mm, m_int)
    return mm .|> Int |> Diagonal
end
