####
#### Init constraints
####

"""
    has_pfinitconstraint(c::ComponentModel)
    has_pfinitconstraint(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has an initialization constraint which depends on the pf state in metadata.

See also: [`get_pfinitconstraints`](@ref), [`set_pfinitconstraint!`](@ref).
"""
has_pfinitconstraint(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfinitconstraint)
has_pfinitconstraint(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfinitconstraint(getcomp(nw, idx))
has_pfinitconstraint(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfinitconstraint(getcomp(nw, idx))

"""
    get_pfinitconstraints(c::NetworkDynamics.ComponentModel)
    get_pfinitconstraints(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the initialization constraints which depend on pf state for the component model.
Returns a vector of constraints, even if only one constraint is present.
May error if no constraints are present. Use `has_pfinitconstraint` to check first.

See also: [`has_pfinitconstraint`](@ref), [`set_pfinitconstraint!`](@ref).
"""
get_pfinitconstraints(c::NetworkDynamics.ComponentModel) = NetworkDynamics._get_list(c, :pfinitconstraint, PFInitConstraint)
get_pfinitconstraints(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfinitconstraints(getcomp(nw, idx))
get_pfinitconstraints(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfinitconstraints(getcomp(nw, idx))

"""
    set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint; check=true)
    set_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint; check=true)

Sets initialization constraints which depend on the powerflow solution to the component.
Accepts a single `PFInitConstraint` or a tuple or vector of them.
Overwrites any existing pf constraints.
See also [`delete_pfinitconstraints!`](@ref), [`add_pfinitconstraint!`](@ref).
"""
function set_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint)
    @nospecialize constraint
    set_metadata!(c, :pfinitconstraint, NetworkDynamics._as_list(PFInitConstraint, constraint, "Constraint"))
end
set_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex, constraint; kw...) = set_pfinitconstraint!(getcomp(nw, idx), constraint; kw...)

"""
    add_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint) -> Bool
    add_pfinitconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint) -> Bool

Adds a new initialization constraint which depends on the powerflow solution to the component.
If constraints already exist, the new constraint is added to the existing ones.
If no constraints exist, this is equivalent to `set_pfinitconstraint!`.

Returns `true` if the constraint was successfully added, `false` if it already exists.

See also [`set_pfinitconstraint!`](@ref), [`delete_pfinitconstraints!`](@ref).
"""
function add_pfinitconstraint!(c::NetworkDynamics.ComponentModel, constraint::PFInitConstraint)
    NetworkDynamics._add_to_list!(c, :pfinitconstraint, PFInitConstraint, constraint)
end
add_pfinitconstraint!(nw::Network, idx::NetworkDynamics.VCIndex, constraint) = add_pfinitconstraint!(getcomp(nw, idx), constraint)
add_pfinitconstraint!(nw::Network, idx::NetworkDynamics.ECIndex, constraint) = add_pfinitconstraint!(getcomp(nw, idx), constraint)

"""
    delete_pfinitconstraints!(c::NetworkDynamics.ComponentModel)
    delete_pfinitconstraints!(nw::Network, idx::Union{VIndex,EIndex})

Removes the powerflow dependent initialization constraint from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the constraint existed and was removed, `false` otherwise.

See also: [`set_pfinitconstraint!`](@ref).
"""
delete_pfinitconstraints!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfinitconstraint)
delete_pfinitconstraints!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfinitconstraints!(getcomp(nw, idx))
delete_pfinitconstraints!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfinitconstraints!(getcomp(nw, idx))

####
#### Init formulas
####

"""
    has_pfinitformula(c::ComponentModel)
    has_pfinitformula(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has an initialization formula which depends on the pf state in metadata.

See also: [`get_pfinitformulas`](@ref), [`set_pfinitformula!`](@ref).
"""
has_pfinitformula(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfinitformula)
has_pfinitformula(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfinitformula(getcomp(nw, idx))
has_pfinitformula(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfinitformula(getcomp(nw, idx))

"""
    get_pfinitformulas(c::NetworkDynamics.ComponentModel)
    get_pfinitformulas(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the initialization formulas which depend on pf state for the component model.
Returns a vector of formulas, even if only one formula is present.
May error if no formulas are present. Use `has_pfinitformula` to check first.

See also: [`has_pfinitformula`](@ref), [`set_pfinitformula!`](@ref).
"""
get_pfinitformulas(c::NetworkDynamics.ComponentModel) = NetworkDynamics._get_list(c, :pfinitformula, PFInitFormula)
get_pfinitformulas(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfinitformulas(getcomp(nw, idx))
get_pfinitformulas(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfinitformulas(getcomp(nw, idx))

"""
    set_pfinitformula!(c::NetworkDynamics.ComponentModel, formula; check=true)
    set_pfinitformula!(nw::Network, idx::Union{VIndex,EIndex}, formula; check=true)

Sets initialization formulas which depend on the powerflow solution to the component.
Accepts a single `PFInitFormula` or a tuple or vector of them.
Overwrites any existing pf formulas.
See also [`delete_pfinitformulas!`](@ref), [`add_pfinitformula!`](@ref).
"""
function set_pfinitformula!(c::NetworkDynamics.ComponentModel, formula)
    @nospecialize formula
    set_metadata!(c, :pfinitformula, NetworkDynamics._as_list(PFInitFormula, formula, "Formula"))
end
set_pfinitformula!(nw::Network, idx::NetworkDynamics.VCIndex, formula; kw...) = set_pfinitformula!(getcomp(nw, idx), formula; kw...)

"""
    add_pfinitformula!(c::NetworkDynamics.ComponentModel, formula::PFInitFormula) -> Bool
    add_pfinitformula!(nw::Network, idx::Union{VIndex,EIndex}, formula) -> Bool

Adds a new initialization formula which depends on the powerflow solution to the component.
If formulas already exist, the new formula is added to the existing ones.
If no formulas exist, this is equivalent to `set_pfinitformula!`.

Returns `true` if the formula was successfully added, `false` if it already exists.

See also [`set_pfinitformula!`](@ref), [`delete_pfinitformulas!`](@ref).
"""
function add_pfinitformula!(c::NetworkDynamics.ComponentModel, formula::PFInitFormula)
    NetworkDynamics._add_to_list!(c, :pfinitformula, PFInitFormula, formula)
end
add_pfinitformula!(nw::Network, idx::NetworkDynamics.VCIndex, formula) = add_pfinitformula!(getcomp(nw, idx), formula)
add_pfinitformula!(nw::Network, idx::NetworkDynamics.ECIndex, formula) = add_pfinitformula!(getcomp(nw, idx), formula)

"""
    delete_pfinitformulas!(c::NetworkDynamics.ComponentModel)
    delete_pfinitformulas!(nw::Network, idx::Union{VIndex,EIndex})

Removes the powerflow dependent initialization formula from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the formula existed and was removed, `false` otherwise.

See also: [`set_pfinitformula!`](@ref).
"""
delete_pfinitformulas!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfinitformula)
delete_pfinitformulas!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfinitformulas!(getcomp(nw, idx))
delete_pfinitformulas!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfinitformulas!(getcomp(nw, idx))

####
#### Powerflow callbacks
####

# Stored as a `Vector{Any}` like the callback lists in NetworkDynamics. A callback is any callable,
# so a tuple or vector is taken as a collection of callbacks and anything else as a single one.
_pfcallback_list(x) = x isa Union{Tuple,AbstractVector} ? collect(Any, x) : Any[x]

"""
    has_pfcallback(c::ComponentModel)
    has_pfcallback(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has powerflow callbacks in metadata.

See also: [`add_pfcallback!`](@ref), [`solve_powerflow`](@ref).
"""
has_pfcallback(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfcallbacks)
has_pfcallback(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfcallback(getcomp(nw, idx))
has_pfcallback(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfcallback(getcomp(nw, idx))

"""
    get_pfcallbacks(c::ComponentModel)
    get_pfcallbacks(nw::Network, idx::Union{VIndex,EIndex})

Gets all powerflow callbacks of the component as a vector, even if there is only a single one.

See also: [`add_pfcallback!`](@ref), [`solve_powerflow`](@ref).
"""
function get_pfcallbacks(c::NetworkDynamics.ComponentModel)
    x = get_metadata(c, :pfcallbacks)
    x isa Vector{Any} ? x : _pfcallback_list(x)
end
get_pfcallbacks(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfcallbacks(getcomp(nw, idx))
get_pfcallbacks(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfcallbacks(getcomp(nw, idx))

"""
    set_pfcallback!(c::ComponentModel, f)
    set_pfcallback!(nw::Network, idx::Union{VIndex,EIndex}, f)

Sets the powerflow callbacks of the component. Overwrites any existing ones.
`f` can be a single callable or a tuple or vector of them.

Powerflow callbacks let iterative powerflow solvers handle things like PV/PQ switching or
stepped shunts. See [`solve_powerflow`](@ref) for details.

See also: [`add_pfcallback!`](@ref), [`delete_pfcallbacks!`](@ref).
"""
function set_pfcallback!(c::NetworkDynamics.ComponentModel, f)
    @nospecialize f
    set_metadata!(c, :pfcallbacks, _pfcallback_list(f))
end
set_pfcallback!(nw::Network, idx::NetworkDynamics.VCIndex, f) = set_pfcallback!(getcomp(nw, idx), f)
set_pfcallback!(nw::Network, idx::NetworkDynamics.ECIndex, f) = set_pfcallback!(getcomp(nw, idx), f)

"""
    add_pfcallback!(c::ComponentModel, f) -> Bool
    add_pfcallback!(nw::Network, idx::Union{VIndex,EIndex}, f) -> Bool

Adds a powerflow callback to the component. Does not overwrite existing ones.
Returns `false` if `f` is already there.

Powerflow callbacks let iterative powerflow solvers handle things like PV/PQ switching or
stepped shunts. See [`solve_powerflow`](@ref) for details.

See also: [`set_pfcallback!`](@ref), [`delete_pfcallbacks!`](@ref).
"""
function add_pfcallback!(c::NetworkDynamics.ComponentModel, f)
    @nospecialize f
    list = has_pfcallback(c) ? get_pfcallbacks(c) : Any[]
    f ∈ list && return false
    set_metadata!(c, :pfcallbacks, push!(list, f))
    return true
end
add_pfcallback!(nw::Network, idx::NetworkDynamics.VCIndex, f) = add_pfcallback!(getcomp(nw, idx), f)
add_pfcallback!(nw::Network, idx::NetworkDynamics.ECIndex, f) = add_pfcallback!(getcomp(nw, idx), f)

"""
    delete_pfcallbacks!(c::ComponentModel)
    delete_pfcallbacks!(nw::Network, idx::Union{VIndex,EIndex})

Removes all powerflow callbacks from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if callbacks existed and were removed, `false` otherwise.

See also: [`set_pfcallback!`](@ref), [`add_pfcallback!`](@ref).
"""
delete_pfcallbacks!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfcallbacks)
delete_pfcallbacks!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfcallbacks!(getcomp(nw, idx))
delete_pfcallbacks!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfcallbacks!(getcomp(nw, idx))
