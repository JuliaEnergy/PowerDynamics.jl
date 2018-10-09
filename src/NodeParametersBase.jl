# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

"Abstract super type for all node parameter types."
abstract type AbstractNodeParameters end
function  Base.show(io::IO, p::AbstractNodeParameters; sep=", ")
    s = "$(nameof(typeof(p)))"
    if !isempty(internalsymbolsof(p))
        s *= "["* join( p |> internalsymbolsof .|> QuoteNode .|> Symbol .|> String, sep) *"]"
    end
    s *= "(" * join(["$(f)=$(repr(getfield(p, f)))" for f in fieldnames(typeof(p))], sep) * ")"
    print(io, s)
end

symbolsof(p::AbstractNodeParameters) = p |> construct_node_dynamics |> symbolsof

abstract type AbstractNodeSymbols end

struct ODENodeSymbols <: AbstractNodeSymbols
    vars::AbstractVector{Symbol}
    dvars::AbstractVector{Symbol}
end
internalsymbolsof(s::ODENodeSymbols) = s.vars
internalsymbolsof(s::AbstractNodeSymbols) = s |> ODENodeSymbols |> internalsymbolsof
internalsymbolsof(s) = s |> symbolsof |> internalsymbolsof
internaldsymbolsof(s::ODENodeSymbols) = s.dvars
internaldsymbolsof(s::AbstractNodeSymbols) = s |> ODENodeSymbols |> internaldsymbolsof
internaldsymbolsof(s) = s |> symbolsof |> internaldsymbolsof
struct DAENodeSymbols <: AbstractNodeSymbols
    odesymbols::ODENodeSymbols
    outvars::AbstractVector{Symbol}
    DAENodeSymbols(s::ODENodeSymbols, outvars) = new(s, outvars)
    DAENodeSymbols(vars, dvars, outvars) = new(ODENodeSymbols(vars, dvars), outvars)
end
ODENodeSymbols(s::DAENodeSymbols) = s.odesymbols
internaloutsymbolsof(s::DAENodeSymbols) = s.outvars
internaloutsymbolsof(s::AbstractNodeSymbols) = s |> DAENodeSymbols |> internaloutsymbolsof
internaloutsymbolsof(s) = s |> symbolsof |> internaloutsymbolsof

Base.convert(::Type{DAENodeSymbols}, s::ODENodeSymbols) = DAENodeSymbols(s, map(x -> Symbol(:out, x), internalsymbolsof(s)))
