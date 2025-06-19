"""
    struct PFInitConstraint{F}
    PFInitConstraint(f, sym, pfsym, dim)

A representation of an additional constraint that is applied during the initialization phase of a component.
In contrast to a [`InitConstraint`](@ref), this constraint may access additonal variables which are available
in the full `NWState` of the solved power flow!

Crucially, this is only necessary for constraints, wnich cannot be expressd in terms of the
**interface variables** (voltages and currents).

See also [`@pfinitconstraint`](@ref) for a macro to create such constraints.
"""
struct PFInitConstraint{F}
    f::F
    sym::Vector{Symbol}
    pfsym::Vector{Symbol}
    dim::Int
    prettyprint::Union{Nothing,String}
end
PFInitConstraint(f, sym, pfsym, dim) = PFInitConstraint(f, sym, pfsym, dim, nothing)
NetworkDynamics.dim(c::PFInitConstraint) = c.dim

# mainly for testing
function (c::PFInitConstraint)(res, u, pfu)
    c.f(res, SymbolicView(u, c.sym), SymbolicView(pfu, c.pfsym))
end

macro pfinitconstraint(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    sym = Symbol[]
    pfsym = Symbol[]
    out = gensym(:out)
    u = gensym(:u)
    pfu = gensym(:pfu)
    body = Expr[]

    dim = 0
    for constraint in ex.args
        constraint isa Union{Expr, QuoteNode} || continue # skip line number nodes
        dim += 1
        wrapped = _wrap_symbols!(constraint, sym, pfsym, u, pfu)
        push!(body, :($(esc(out))[$dim] = $wrapped))
    end
    unique!(sym)
    unique!(pfsym)

    s = join(string.(body), "\n")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(out))\"))))" => "    out")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(u))\"))))" => "u")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(pfu))\"))))" => "pfu")
    s = "PFInitConstraint($sym, $pfsym, $dim) do out, u, pfu\n" * s * "\nend"

    quote
        PFInitConstraint($sym, $pfsym, $dim, $s) do $(esc(out)), $(esc(u)), $(esc(pfu))
            $(body...)
            nothing
        end
    end
end
function _wrap_symbols!(ex, sym, pfsym, u, pfu)
    postwalk(ex) do x
        if x isa QuoteNode && x.value isa Symbol
            push!(sym, x.value)
            :($(esc(u))[$x])
        elseif @capture(x, @pf($(esc(u))[sc_]))
            @assert sc isa QuoteNode && sc.value isa Symbol
            # delete sc from symbols
            deleteat!(sym, findlast(isequal(sc.value), sym))
            push!(pfsym, sc.value)
            :($(esc(pfu))[$sc])
        else
            x
        end
    end
end

function pfinitc_to_initc(pfic::PFInitConstraint, pfstate::NWState, cidx)
    VEIndex = NetworkDynamics._baseT(cidx)
    @assert cidx isa VEIndex{<:Any, Nothing}
    pfvec = pfstate[collect(VEIndex(cidx.compidx, pfic.pfsym))]

    pfu = SymbolicView(pfvec, pfic.sym)
    f = pfic.f

    InitConstraint(pfic.sym, pfic.dim) do res, u
        f(res, u, pfu)
    end
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::PFInitConstraint))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
    end
end
