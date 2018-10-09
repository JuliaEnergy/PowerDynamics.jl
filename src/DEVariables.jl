# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: length, view, @__doc__
using MacroTools

issymbol(::Symbol) = true
issymbol(::Any) = false

"Abstract super type for all variables that AbstractNodeDynamics sub types can be called with."
abstract type AbstractDEVariable end

"""
    function mapfields(f, s, args...)

Applies `f` to all fields of (the struct) `s` giving `args...` as additional arguments.
"""
@generated function mapfields(f, s, args...)
    Expr(:call, nameof(s), [:(f(s.$field, args...)) for field in fieldnames(s) ]...)
end

"Extend view from arrays to subtypes of [`AbstractDEVariable`](#ref)."
function view(var::AbstractDEVariable, range)
    mapfields(view, var, range)
end

"Extend complexview from arrays to subtypes of [`AbstractDEVariable`](#ref)."
function complexview(var::AbstractDEVariable, i0, num)
    mapfields(complexview, var, i0, num)
end

length(var::AbstractDEVariable) = length(var.val)

_excomparison(ex) = [ex]
_excomparison(ex, exs...) = [[ex, :(==)]; _excomparison(exs...)]
"Create an expresseion where `==` is applied between all the expressions given as argument here."
excomparison(exs) = Expr(:comparison, _excomparison(exs...)...)

"""
Basically, this macro generates all the constructors (internal and external) for a subtype of AbstracDEVariable.

If you really want to understand this macro, uncomment the `println(ex)` statement at the end and check the resulting
generated expression.
"""
function create_DEVariable(structdef, outputfield::Symbol)
    structdef = deepcopy(structdef)
    @assert MacroTools.isstructdef(structdef)
    @assert @capture(structdef, struct name_{elementtypevars__} <: supername_
        typedfields__
    end)
    @capture(structdef, struct typedname_ <: _
        __
    end)
    @assert length(typedfields) >= 1
    fields = [field.args[1] for field in typedfields]
    typedinputfields = [field for field in typedfields if field.args[1] != outputfield]
    inputfields = [field.args[1] for field in typedinputfields]
    basic_constructor = Expr(:function, Expr(:where, Expr(:call, typedname, typedfields...), elementtypevars...), quote
            $( length(fields) > 1 ? :(@assert $(excomparison(map(field -> Expr(:call, :length, field), fields)))) : nothing)
            $(Expr(:call, :new, fields...))
        end
    )
    external_basic_constructor = :($(Expr(:where, Expr(:call, name, typedfields...), elementtypevars...)) = $(Expr(:call, typedname, fields...)))
    external_add_constructor = :($(Expr(:where, Expr(:call, name, typedinputfields...), elementtypevars...)) = $(Expr(:call, name, typedinputfields..., :(similar($(inputfields[1]))) )) )
    kw_constructor = Expr(:function, Expr(:call, name, Expr(:parameters, inputfields..., Expr(:kw, outputfield, :(nothing)))), quote
            if $outputfield == nothing
                $(Expr(:call, name, inputfields...))
            else
                $(Expr(:call, name, fields...))
            end
        end
    )
    append!(structdef.args[3].args, [basic_constructor])
    ex = quote
        @__doc__ $structdef
        $external_basic_constructor
        $external_add_constructor
        $kw_constructor
    end
    ex
end

macro DEVariable(structdef, outputfield)
    mainex = create_DEVariable(structdef, outputfield)
    @capture(structdef, struct name_{__} <: _
        __
    end)
    mainexstr = "$(copy(mainex)|>rlr)"
    showex = :(showdefinition(io::IO, ::Type{$name}) = println(io, $mainexstr))
    append!(mainex.args, [showex])
    return esc(mainex)
end

"Abstract super type for all Variables for ODE-type node dynamics."
abstract type AbstractODEVariable <: AbstractDEVariable end
"Variables for ODE-type node dynamics."
@DEVariable struct ODEVariable{Tval, Tddt} <: AbstractODEVariable
    val::AbstractVector{Tval}
    ddt::AbstractVector{Tddt}
end ddt

"Abstract super type for all Variables for DAE-type node dynamics."
abstract type AbstractDAEVariable <: AbstractDEVariable end
"Variables for DAE-type node dynamics."
@DEVariable struct DAEVariable{Tval, Tddt, Tout} <: AbstractDAEVariable
    val::AbstractVector{Tval}
    ddt::AbstractVector{Tddt}
    out::AbstractVector{Tout}
end out
