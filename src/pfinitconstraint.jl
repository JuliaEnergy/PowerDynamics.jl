"""
    struct PFInitConstraint{F}
    PFInitConstraint(f, sym, pfsym, dim)

A representation of an additional constraint that is applied during the initialization phase of a component.
In contrast to a [`InitConstraint`](@extref NetworkDynamics.InitConstraint), this constraint may access additional variables which are available
in the full `NWState` of the solved power flow!

Crucially, this is only necessary for constraints, which cannot be expressed in terms of the
**interface variables** (voltages and currents).

See also: [`@pfinitconstraint`](@ref) for a macro to create such constraints,
[`PFInitConstraint`](@ref), [`set_pfinitconstraint!`](@ref),
[`add_pfinitconstraint!`](@ref)
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

function Base.:(==)(a::PFInitConstraint, b::PFInitConstraint)
    typeof(a) == typeof(b) && NetworkDynamics.equal_fields(a, b)
end

# mainly for testing
function (c::PFInitConstraint)(res, u, pfu)
    c.f(res, SymbolicView(u, c.sym), SymbolicView(pfu, c.pfsym))
end

"""
    struct PFInitFormula{F}
    PFInitFormula(f, outsym, sym, pfsym)

A representation of an initialization formula that is applied during the initialization phase of a component.
In contrast to a [`InitFormula`](@extref NetworkDynamics.InitFormula), this formula may access additional variables which are available
in the full `NWState` of the solved power flow!

Crucially, this is only necessary for formulas, which cannot be expressed in terms of the
**interface variables** (voltages and currents).

Similar to InitFormula, this sets defaults rather than adding constraint equations.
The formula is applied early in the initialization pipeline before constraints are solved.

See also: [`@pfinitformula`](@ref) for a macro to create such formulas,
[`PFInitFormula`](@ref), [`set_pfinitformula!`](@ref),
[`add_pfinitformula!`](@ref)
"""
struct PFInitFormula{F}
    f::F
    outsym::Vector{Symbol}
    sym::Vector{Symbol}
    pfsym::Vector{Symbol}
    prettyprint::Union{Nothing,String}
end
PFInitFormula(f, outsym, sym, pfsym) = PFInitFormula(f, outsym, sym, pfsym, nothing)

function Base.:(==)(a::PFInitFormula, b::PFInitFormula)
    typeof(a) == typeof(b) && NetworkDynamics.equal_fields(a, b)
end

# mainly for testing
function (c::PFInitFormula)(res, u, pfu)
    c.f(SymbolicView(res, c.outsym), SymbolicView(u, c.sym), SymbolicView(pfu, c.pfsym))
end

"""
    @pfinitconstraint expr
    @pfinitconstraint begin
        constraint1
        constraint2
    end

Create a [`PFInitConstraint`](@ref) using macro syntax. Component variables are accessed with
`:symbol` and power flow state variables with `@pf(:symbol)`. Multiple constraints can be
defined in a begin...end block.

For example, the following code creates constraints `Pset`/`Qset` of the dynamic model
in reference to the variables `P` and `Q` of the powerflow model.

```julia
@pfinitconstraint :dynamic₊Pset - @pf(:static₊P)
@pfinitconstraint begin
    :dynamic₊Pset - @pf(:static₊P)
    :dynamic₊Qset - @pf(:static₊Q)
end
```

Note that this particiular example could have been a `PFInitFormula`, since we can solve it directly.

See also: [`PFInitConstraint`](@ref), [`set_pfinitconstraint!`](@ref), [`add_pfinitconstraint!`](@ref)
"""
macro pfinitconstraint(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    # capture the macro-form source before the symbols get wrapped into `u[:sym]`/`pfu[:sym]`
    s = _macro_source_string("@pfinitconstraint", ex)

    sym = Symbol[]
    pfsym = Symbol[]
    out = :__out__
    u = :__u__
    pfu = :__pfu__
    body = Expr[]

    dim = 0
    for constraint in ex.args
        constraint isa Union{Expr, QuoteNode} || continue # skip line number nodes
        dim += 1
        wrapped = _wrap_symbols!(constraint, sym, pfsym, u, pfu)
        push!(body, :($out[$dim] = $wrapped))
    end
    unique!(sym)
    unique!(pfsym)

    body_esc = _escape_all.(body)

    quote
        PFInitConstraint($sym, $pfsym, $dim, $s) do $(esc(out)), $(esc(u)), $(esc(pfu))
            $(body_esc...)
            nothing
        end
    end
end
function _wrap_symbols!(ex, sym, pfsym, u, pfu)
    postwalk(ex) do x
        if x isa QuoteNode && x.value isa Symbol
            push!(sym, x.value)
            :($u[$x])
        elseif @capture(x, @pf($u[sc_]))
            @assert sc isa QuoteNode && sc.value isa Symbol
            # delete sc from symbols (it was added by the QuoteNode branch above)
            deleteat!(sym, findlast(isequal(sc.value), sym))
            push!(pfsym, sc.value)
            :($pfu[$sc])
        else
            x
        end
    end
end
function _escape_all(ex::Expr)
    postwalk(ex) do x
        x isa Symbol ? esc(x) : x
    end
end

# Reproduce the `@name begin … end` source the user wrote, so a constraint/formula prints as
# its macro call rather than as the expanded inline `do`-block — it reads nicer. Must be called
# before the symbols are wrapped into `u[:sym]`/`pfu[:sym]`. The `#= file:line =#` location
# comment on nested macro calls (`@pf(:x)`) survives `remove_linenums!`, so strip it here.
function _macro_source_string(macroname, ex)
    lines = [replace(string("    ", a), r"#=.*?=# " => "") for a in ex.args if a isa Union{Expr, QuoteNode}]
    string(macroname, " begin\n", join(lines, "\n"), "\nend")
end

"""
    @pfinitformula expr
    @pfinitformula begin
        :var1 = expr1
        :var2 = expr2
    end

Create a [`PFInitFormula`](@ref) using macro syntax. Component variables are accessed with
`:symbol` and power flow state variables with `@pf(:symbol)`. Multiple formulas can be
defined in a begin...end block.

Unlike constraints, formulas use assignment syntax (`:var = expression`) to set variable values
during initialization. The left-hand side specifies output variables, and the right-hand side
can access both component variables and power flow state variables.

For example, the following code creates formulas for `Pset`/`Qset` of the dynamic model
in reference to the variables `P` and `Q` of the powerflow model.

```julia
@pfinitformula :dynamic₊Pset = @pf(:static₊P)
@pfinitformula begin
    :dynamic₊Pset = @pf(:static₊P)
    :dynamic₊Qset = @pf(:static₊Q)
end
```

See also: [`PFInitFormula`](@ref), [`set_pfinitformula!`](@ref), [`add_pfinitformula!`](@ref)
"""
macro pfinitformula(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    # capture the macro-form source before the symbols get wrapped into `u[:sym]`/`pfu[:sym]`
    s = _macro_source_string("@pfinitformula", ex)

    sym = Symbol[]
    pfsym = Symbol[]
    outsym = Symbol[]
    out = :__out__
    u = :__u__
    pfu = :__pfu__
    body = Expr[]

    for formula in ex.args
        formula isa Union{Expr, QuoteNode} || continue # skip line number nodes

        # Parse assignment expressions
        if formula isa Expr && formula.head == :(=)
            lhs = formula.args[1]
            rhs = formula.args[2]

            # Extract output symbol from left-hand side
            if lhs isa QuoteNode && lhs.value isa Symbol
                push!(outsym, lhs.value)
                wrapped_rhs = _wrap_symbols!(rhs, sym, pfsym, u, pfu)
                push!(body, :($out[$lhs] = $wrapped_rhs))
            else
                error("Left-hand side of formula assignment must be a quoted symbol like :var")
            end
        else
            error("Formula expressions must be assignments of the form :var = expression")
        end
    end

    unique!(sym)
    unique!(pfsym)
    unique!(outsym)

    body_esc = _escape_all.(body)

    quote
        PFInitFormula($outsym, $sym, $pfsym, $s) do $(esc(out)), $(esc(u)), $(esc(pfu))
            $(body_esc...)
            nothing
        end
    end
end

"""
    specialize_pfinitconstraints(nw, pfs)

Convert all [`PFInitConstraint`](@ref)s in the network to regular [`InitConstraint`](@ref)s
by specializing them with the power flow solution. Scans through all components and extracts
power flow variables needed by each constraint, creating specialized versions with those
values embedded.

Called internally by [`initialize_from_pf[!]`](@ref initialize_from_pf).
"""
function specialize_pfinitconstraints(nw, pfs)
    dict = Dict{NetworkDynamics.SymbolicIndex, Vector{InitConstraint}}()
    vidxs = (VIndex(i) for i in 1:nv(nw))
    eidxs = (EIndex(i) for i in 1:ne(nw))
    for cidx in Iterators.flatten((vidxs, eidxs))
        c = nw[cidx]
        if has_pfinitconstraint(c)
            pfics = get_pfinitconstraints(c)
            ics = [specialize_pfinitconstraint(pfic, pfs, cidx) for pfic in pfics]
            dict[cidx] = ics
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
    VEIndex = NetworkDynamics.idxtype(cidx)
    @assert cidx isa VEIndex{<:Any, Nothing}
    pfvec = pfstate[collect(VEIndex(cidx.compidx, pfic.pfsym))]

    pfu = SymbolicView(pfvec, pfic.pfsym)
    f = pfic.f

    InitConstraint(pfic.sym, pfic.dim) do res, u
        f(res, u, pfu)
    end
end

"""
    specialize_pfinitformulas(nw, pfs)

Convert all [`PFInitFormula`](@ref)s in the network to regular [`InitFormula`](@ref)s
by specializing them with the power flow solution. Scans through all components and extracts
power flow variables needed by each formula, creating specialized versions with those
values embedded.

Called internally by [`initialize_from_pf[!]`](@ref initialize_from_pf).
"""
function specialize_pfinitformulas(nw, pfs)
    dict = Dict{NetworkDynamics.SymbolicIndex, Vector{InitFormula}}()
    vidxs = (VIndex(i) for i in 1:nv(nw))
    eidxs = (EIndex(i) for i in 1:ne(nw))
    for cidx in Iterators.flatten((vidxs, eidxs))
        c = nw[cidx]
        if has_pfinitformula(c)
            pfifs = get_pfinitformulas(c)
            ifs = [specialize_pfinitformula(pfif, pfs, cidx) for pfif in pfifs]
            dict[cidx] = ifs
        end
    end
    dict
end

"""
    specialize_pfinitformula(pfif::PFInitFormula, pfstate::NWState, cidx)

Convert a single [`PFInitFormula`](@ref) to a regular [`InitFormula`](@ref) by extracting
the required power flow variables from `pfstate` and embedding them in the formula function.

Called by [`specialize_pfinitformulas`](@ref) for each component with a PFInitFormula.
"""
function specialize_pfinitformula(pfif::PFInitFormula, pfstate::NWState, cidx)
    VEIndex = NetworkDynamics.idxtype(cidx)
    @assert cidx isa VEIndex{<:Any, Nothing}
    pfvec = pfstate[collect(VEIndex(cidx.compidx, pfif.pfsym))]

    pfu = SymbolicView(pfvec, pfif.pfsym)
    f = pfif.f

    InitFormula(pfif.outsym, pfif.sym) do res, u
        f(res, u, pfu)
    end
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::PFInitFormula))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
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

"""
    copy_pf_parameters(cm::ComponentModel) -> PFInitFormula

Creates a [`PFInitFormula`](@ref) that copies all parameters from the powerflow model
to the component model. This formula can then be added to the component using
[`add_pfinitformula!`](@ref).

This is useful for components where the powerflow and dynamic models should have
identical parameter values, ensuring consistency between the two models.

See also: [`PFInitFormula`](@ref), [`add_pfinitformula!`](@ref)
"""
function copy_pf_parameters(cm::NetworkDynamics.ComponentModel)
    cmpf = powerflow_model(cm)
    if Set(psym(cm)) != Set(psym(cmpf))
        throw(ArgumentError("Cannot add parameter copy from powerflow model to component :$(cm.name)! The powerflow model has different parameters than the component model!"))
    end
    let _sym=psym(cm), _pfsym=psym(cmpf)
        PFInitFormula(_sym, Symbol[], _pfsym) do res, _, pfu
            for p in _sym
                res[p] = pfu[p]
            end
            nothing
        end
    end
end
