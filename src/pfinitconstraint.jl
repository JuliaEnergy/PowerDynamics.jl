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

"""
    @pfinitconstraint expr
    @pfinitconstraint begin
        constraint1
        constraint2
    end

Create a [`PFInitConstraint`](@ref) using macro syntax. Component variables are accessed with
`:symbol` and power flow state variables with `@pf :symbol`. Multiple constraints can be
defined in a begin...end block.

See also: [`PFInitConstraint`](@ref), [`set_pfinitconstraint!`](@ref)
"""
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

"""
    specialize_pfinitconstraints(nw, pfs)

Convert all [`PFInitConstraint`](@ref)s in the network to regular [`InitConstraint`](@ref)s
by specializing them with the power flow solution. Scans through all components and extracts
power flow variables needed by each constraint, creating specialized versions with those
values embedded.

Called internally by [`initialize_from_pf[!]`](@ref).
"""
function specialize_pfinitconstraints(nw, pfs)
    dict = Dict{NetworkDynamics.SymbolicIndex, InitConstraint}()
    vidxs = (VIndex(i) for i in 1:nv(nw))
    eidxs = (EIndex(i) for i in 1:ne(nw))
    for cidx in Iterators.flatten((vidxs, eidxs))
        c = nw[cidx]
        if has_pfinitconstraint(c)
            pfic = get_pfinitconstraint(c)
            ic = specialize_pfinitconstraint(pfic, pfs, cidx)
            dict[cidx] = pfic
        end
    end
    dict
end

"""
    specialize_pfinitconstraint(pfic::PFInitConstraint, pfstate::NWState, cidx)

Convert a single [`PFInitConstraint`](@ref) to a regular [`InitConstraint`](@ref) by extracting
the required power flow variables from `pfstate` and embedding them in the constraint function.

Called by [`specialize_pfinitconstraints`](@ref) for each component with a PFInitConstraint.
"""
function specialize_pfinitconstraint(pfic::PFInitConstraint, pfstate::NWState, cidx)
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

####
#### Init constraints
####

"""
    has_pfinitconstraint(c::ComponentModel)
    has_pfinitconstraint(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has an initialization constraint which depends on the pf state in metadata.

See also: [`get_pfinitconstraint`](@ref), [`set_pfinitconstraint!`](@ref).
"""
has_pfinitconstraint(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfinitconstraint)
has_pfinitconstraint(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfinitconstraint(getcomp(nw, idx))
has_pfinitconstraint(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfinitconstraint(getcomp(nw, idx))

"""
    get_pfinitconstraint(c::NetworkDynamics.ComponentModel)
    get_pfinitconstraint(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the initialization constraint which depends on pf state for the component model.
May error if no constraint is present. Use `has_pfinitconstraint` to check first.

See also: [`has_pfinitconstraint`](@ref), [`set_pfinitconstraint!`](@ref).
"""
get_pfinitconstraint(c::NetworkDynamics.ComponentModel) = get_metadata(c, :pfinitconstraint)::PFInitConstraint
get_pfinitconstraint(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfinitconstraint(getcomp(nw, idx))
get_pfinitconstraint(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfinitconstraint(getcomp(nw, idx))

"""
    set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint; check=true)
    set_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint; check=true)

Sets an additional initialization constraint which depends on the powerflow solution to
the component. Overwrites any existing pf constraints.
See also [`delete_pfinitconstraint!`](@ref).
"""
function set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint)
    set_metadata!(c, :pfinitconstraint, constraint)
end
set_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex, constraint; kw...) = set_pfinitconstraint(getcomp(nw, idx), constraint; kw...)

"""
    delete_pfinitconstraint!(c::NetworkDynamics.ComponentModel)
    delete_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex})

Removes the powerflow dependent initialization constraint from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the constraint existed and was removed, `false` otherwise.

See also: [`set_pfinitconstraint!`](@ref).
"""
delete_pfinitconstraint!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfinitconstraint)
delete_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfinitconstraint!(getcomp(nw, idx))
delete_pfinitconstraint!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfinitconstraint!(getcomp(nw, idx))
