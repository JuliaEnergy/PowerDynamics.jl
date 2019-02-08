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
