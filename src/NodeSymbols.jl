
"""to be documented"""
function symbolsof end

symbolsof(::N) where {N <: AbstractNodeParameters} = symbolsof(N)
symbolsof(::AbstractNodeDynamics{N}) where {N} = symbolsof(N)

abstract type AbstractNodeSymbols end

struct ODENodeSymbols <: AbstractNodeSymbols
    vars::AbstractVector{Symbol}
    dvars::AbstractVector{Symbol}
end
"Get the symbols representing the internal variables of the node."
internalsymbolsof(s::ODENodeSymbols) = s.vars
internalsymbolsof(s::AbstractNodeSymbols) = s |> ODENodeSymbols |> internalsymbolsof
internalsymbolsof(s) = s |> symbolsof |> internalsymbolsof
"Get the symbols representing the derivative of the internal variables of the node."
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
"Get the symbols representing the output of the internal variables of the node."
internaloutsymbolsof(s::DAENodeSymbols) = s.outvars
internaloutsymbolsof(s::AbstractNodeSymbols) = [Symbol(:out, sym) for sym in internalsymbolsof(s)]
internaloutsymbolsof(s) = s |> symbolsof |> internaloutsymbolsof

Base.convert(::Type{DAENodeSymbols}, s::ODENodeSymbols) = DAENodeSymbols(s, map(x -> Symbol(:out, x), internalsymbolsof(s)))
