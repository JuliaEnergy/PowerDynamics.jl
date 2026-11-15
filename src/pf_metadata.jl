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
