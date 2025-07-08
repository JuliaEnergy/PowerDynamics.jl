"""
    struct PFInitConstraint{F}
    PFInitConstraint(f, sym, pfsym, dim)

A representation of an additional constraint that is applied during the initialization phase of a component.
In contrast to a [`InitConstraint`](@ref), this constraint may access additional variables which are available
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
In contrast to a [`InitFormula`](@ref), this formula may access additional variables which are available
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
`:symbol` and power flow state variables with `@pf :symbol`. Multiple constraints can be
defined in a begin...end block.

See also: [`PFInitConstraint`](@ref), [`set_pfinitconstraint!`](@ref), [`add_pfinitconstraint!`](@ref)
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
    @pfinitformula expr
    @pfinitformula begin
        :var1 = expr1
        :var2 = expr2
    end

Create a [`PFInitFormula`](@ref) using macro syntax. Component variables are accessed with
`:symbol` and power flow state variables with `@pf :symbol`. Multiple formulas can be
defined in a begin...end block.

Unlike constraints, formulas use assignment syntax (`:var = expression`) to set variable values
during initialization. The left-hand side specifies output variables, and the right-hand side
can access both component variables and power flow state variables.

See also: [`PFInitFormula`](@ref), [`set_pfinitformula!`](@ref), [`add_pfinitformula!`](@ref)
"""
macro pfinitformula(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    sym = Symbol[]
    pfsym = Symbol[]
    outsym = Symbol[]
    out = gensym(:out)
    u = gensym(:u)
    pfu = gensym(:pfu)
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
                push!(body, :($(esc(out))[$lhs] = $wrapped_rhs))
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

    s = join(string.(body), "\n")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(out))\"))))" => "    out")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(u))\"))))" => "u")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(pfu))\"))))" => "pfu")
    s = "PFInitFormula($outsym, $sym, $pfsym) do out, u, pfu\n" * s * "\nend"

    quote
        PFInitFormula($outsym, $sym, $pfsym, $s) do $(esc(out)), $(esc(u)), $(esc(pfu))
            $(body...)
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

Called internally by [`initialize_from_pf[!]`](@ref).
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
    VEIndex = NetworkDynamics._baseT(cidx)
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

Called internally by [`initialize_from_pf[!]`](@ref).
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
    VEIndex = NetworkDynamics._baseT(cidx)
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
    get_pfinitconstraints(c::NetworkDynamics.ComponentModel)
    get_pfinitconstraints(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the initialization constraints which depend on pf state for the component model.
Returns a tuple of constraints, even if only one constraint is present.
May error if no constraints are present. Use `has_pfinitconstraint` to check first.

See also: [`has_pfinitconstraint`](@ref), [`set_pfinitconstraint!`](@ref).
"""
function get_pfinitconstraints(c::NetworkDynamics.ComponentModel)
    constraint = get_metadata(c, :pfinitconstraint)
    return constraint isa Tuple ? constraint : (constraint,)
end
get_pfinitconstraints(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfinitconstraints(getcomp(nw, idx))
get_pfinitconstraints(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfinitconstraints(getcomp(nw, idx))

"""
    set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint; check=true)
    set_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint; check=true)

Sets initialization constraints which depend on the powerflow solution to the component.
Accepts either a single `PFInitConstraint` or a tuple of `PFInitConstraint` objects.
Overwrites any existing pf constraints.
See also [`delete_pfinitconstraint!`](@ref), [`add_pfinitconstraint!`](@ref).
"""
function set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::Union{PFInitConstraint, Tuple{Vararg{PFInitConstraint}}})
    set_metadata!(c, :pfinitconstraint, constraint)
end
set_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex, constraint; kw...) = set_pfinitconstraint!(getcomp(nw, idx), constraint; kw...)

"""
    add_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint) -> Bool
    add_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint) -> Bool

Adds a new initialization constraint which depends on the powerflow solution to the component.
If constraints already exist, the new constraint is added to the existing ones.
If no constraints exist, this is equivalent to `set_pfinitconstraint!`.

Returns `true` if the constraint was successfully added, `false` if it already exists.

See also [`set_pfinitconstraint!`](@ref), [`delete_pfinitconstraint!`](@ref).
"""
function add_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint)
    if has_pfinitconstraint(c)
        existing_constraints = get_pfinitconstraints(c)

        constraint ∈ existing_constraints && return false

        new_constraints = (existing_constraints..., constraint)
        set_metadata!(c, :pfinitconstraint, new_constraints)
    else
        set_metadata!(c, :pfinitconstraint, constraint)
    end
    return true
end
add_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex, constraint) = add_pfinitconstraint!(getcomp(nw, idx), constraint)
add_pfinitconstraint!(nw::Network, idx::NetworkDynamics.ECIndex, constraint) = add_pfinitconstraint!(getcomp(nw, idx), constraint)

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

####
#### Init formulas
####

"""
    has_pfinitformula(c::ComponentModel)
    has_pfinitformula(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has an initialization formula which depends on the pf state in metadata.

See also: [`get_pfinitformula`](@ref), [`set_pfinitformula!`](@ref).
"""
has_pfinitformula(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfinitformula)
has_pfinitformula(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfinitformula(getcomp(nw, idx))
has_pfinitformula(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfinitformula(getcomp(nw, idx))

"""
    get_pfinitformulas(c::NetworkDynamics.ComponentModel)
    get_pfinitformulas(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the initialization formulas which depend on pf state for the component model.
Returns a tuple of formulas, even if only one formula is present.
May error if no formulas are present. Use `has_pfinitformula` to check first.

See also: [`has_pfinitformula`](@ref), [`set_pfinitformula!`](@ref).
"""
function get_pfinitformulas(c::NetworkDynamics.ComponentModel)
    formula = get_metadata(c, :pfinitformula)
    return formula isa Tuple ? formula : (formula,)
end
get_pfinitformulas(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfinitformulas(getcomp(nw, idx))
get_pfinitformulas(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfinitformulas(getcomp(nw, idx))

"""
    set_pfinitformula!(c::NetworkDynamics.ComponentModel, formula; check=true)
    set_pfinitformula!(nw::Network, idx::Union{VIndex,EIndex}, formula; check=true)

Sets initialization formulas which depend on the powerflow solution to the component.
Accepts either a single `PFInitFormula` or a tuple of `PFInitFormula` objects.
Overwrites any existing pf formulas.
See also [`delete_pfinitformula!`](@ref), [`add_pfinitformula!`](@ref).
"""
function set_pfinitformula!(c::NetworkDynamics.ComponentModel, formula::Union{PFInitFormula, Tuple{Vararg{PFInitFormula}}})
    set_metadata!(c, :pfinitformula, formula)
end
set_pfinitformula!(nw::Network, idx::NetworkDynamics.VCIndex, formula; kw...) = set_pfinitformula!(getcomp(nw, idx), formula; kw...)

"""
    add_pfinitformula!(c::NetworkDynamics.ComponentModel, formula::PFInitFormula) -> Bool
    add_pfinitformula!(nw::Network, idx::Union{VIndex,EIndex}, formula) -> Bool

Adds a new initialization formula which depends on the powerflow solution to the component.
If formulas already exist, the new formula is added to the existing ones.
If no formulas exist, this is equivalent to `set_pfinitformula!`.

Returns `true` if the formula was successfully added, `false` if it already exists.

See also [`set_pfinitformula!`](@ref), [`delete_pfinitformula!`](@ref).
"""
function add_pfinitformula!(c::NetworkDynamics.ComponentModel, formula::PFInitFormula)
    if has_pfinitformula(c)
        existing_formulas = get_pfinitformulas(c)

        formula ∈ existing_formulas && return false

        new_formulas = (existing_formulas..., formula)
        set_metadata!(c, :pfinitformula, new_formulas)
    else
        set_metadata!(c, :pfinitformula, formula)
    end
    return true
end
add_pfinitformula!(nw::Network, idx::NetworkDynamics.VCIndex, formula) = add_pfinitformula!(getcomp(nw, idx), formula)
add_pfinitformula!(nw::Network, idx::NetworkDynamics.ECIndex, formula) = add_pfinitformula!(getcomp(nw, idx), formula)

"""
    delete_pfinitformula!(c::NetworkDynamics.ComponentModel)
    delete_pfinitformula!(nw::Network, idx::Union{VIndex,EIndex})

Removes the powerflow dependent initialization formula from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the formula existed and was removed, `false` otherwise.

See also: [`set_pfinitformula!`](@ref).
"""
delete_pfinitformula!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfinitformula)
delete_pfinitformula!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfinitformula!(getcomp(nw, idx))
delete_pfinitformula!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfinitformula!(getcomp(nw, idx))

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
