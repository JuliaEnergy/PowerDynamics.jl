# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: convert

# G = grid type
# V = Value Type of entries (e.g. Float64)
# T = type of the time entry
abstract type AbstractState{G<: GridDynamics, V, T} end

@doc doc"""
```Julia

    BaseState(grid, vec)

```

Encode a state vector and the corresponding rhs information.

# Keyword Arguments
- `grid` is a [`GridDynamics`](@ref) instance that contains the overall system rhs.
- `vec` is a state vector of the system who's length is given by the total
        number of internal and voltage variables.

# Indexing

In an instance `b` of of a `BaseState` behaves like an
[`Array`](@ref), i.e. you can access the ``j``-th element
of the state vector (and set it to a value ``ξ``) by calling `b[j] ( = ξ )`.

"""
struct BaseState{G <: GridDynamics, V} #<: AbstractVector{V} where {G <: GridDynamics}
    grid::G
    vec::AbstractVector{V}
    function BaseState{G, V}(grid::G, vec::AbstractVector{V}) where {G <: GridDynamics, V}
        @assert grid.rhs.systemsize == length(vec)
        new{G, V}(grid, vec)
    end
end
BaseState(grid::G, vec::AbstractVector{V}) where {G <: GridDynamics, V} = BaseState{G, V}(grid, vec)

Base.size(s::BaseState) = size(s.vec)
Base.getindex(s::BaseState, n) = getindex(s.vec, n)
Base.setindex!(s::BaseState, v,  n) = setindex!(s.vec, v, n)
Base.copy(s::BaseState) = BaseState(GridDynamics(s), deepcopy(s.vec)) # grid dynamics should not be copied
GridDynamics(s::BaseState) = s.grid

convert(::Type{<:AbstractVector{V}}, s::BaseState{G, V}) where {G, V} = s.vec

@doc doc"""
```Julia

    State(base, t)
    State(grid, vec, t)

```

Encode the information on the value of a state vector at a particular time point.

# Keyword Arguments
- `base` is an instance of a `BaseState`, essentially it contains the state
    vector and the complete rhs of the system. Instead of `base`, you can also
    directly use a [`GridDynamics`](@ref) instance `grid` and a properly sized
    state vector `vec` to instantiate a `State`.
- `t` is a time point associated to the `base`. It defaults to `nothing`.

# Indexing

Concerning the indexing, a `State` object ``s`` basically behaves like a
an [`Array`](@ref).
There are plenty of convenient ways to access its contents at a node ``j``
by using a particular [`Symbol`](@ref):

* `s[j, :u]`: complex voltage
* `s[j, :v]`: voltage magnitude
* `s[j, :φ]`: voltage angle
* `s[j, :i]`: complex nodal current
* `s[j, :iabs]`: nodal current magnitude
* `s[j, :δ]`: nodal current angle
* `s[j, :s]`: apparent power
* `s[j, :p]`: real power
* `s[j, :q]`: reactive power

Currently, setting the state value is only implemented for ``u`` and ``v``,
the other quantities are derived automatically.

When a node ``j`` has internal variables, you can access (and set) the ``k``-th internal
variable by calling

`s[j, :int, k]`.

The internal variables can be also directly accessed with symbols, i.e.

`s[j, :ω]`

returns the frequency ``ω`` at node ``j``.
To find out the proper symbol, the easiest way is to look into the docs
of the corresponding [`AbstractNodeDynamics`](@ref) subtype or simply at the
output of `print`:

    julia> println(SwingEq(H=2, P=3, D=4, Ω=5))
    SwingEq[:ω](H=2, P=3, D=4, Ω=5)

"""
struct State{G, V, T} <: AbstractState{G, V, T}
    base::BaseState{G, V}
    t::T # keep the time here, so it's easier to define time series later
end
State(b::BaseState; t=nothing) = State(b, t)
State(g::GridDynamics, v::AbstractVector; t=nothing) = State(BaseState(g, v), t)

# Interface Functions
BaseState(s::State) = s.base
Base.copy(s::State) = State(copy(s.base), deepcopy(s.t)) # grid dynamics should not be copied


# derived function for all AbstractState
# Base.size(s::AbstractState) = size(BaseState(s))
# Base.size(s::AbstractState) = (length(Nodes(s)), Set([:u, :v, :]))
GridDynamics(s::AbstractState) = GridDynamics(BaseState(s))
Base.getindex(s::AbstractState, n::Colon, sym::Symbol, args...) = getindex(s, eachindex(Nodes(s)), sym, args...)
Base.getindex(s::AbstractState, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(Nodes(s))  )
        throw(BoundsError(s, n))
    end
    getindex(s, n, Val{sym}, args...)
end
Base.getindex(s::AbstractState, n, ::Type{Val{:u}}) = getindex(BaseState(s), 2 .* n .- 1) + im .* getindex(BaseState(s), 2 .* n)
Base.getindex(s::AbstractState, n, ::Type{Val{:v}}) = abs.(s[n, :u])
Base.getindex(s::AbstractState, n, ::Type{Val{:φ}}) = angle.(s[n, :u])
Base.getindex(s::AbstractState, n, ::Type{Val{:i}}) = getindex(AdmittanceLaplacian(s) * s[:, :u], n)
Base.getindex(s::AbstractState, n, ::Type{Val{:iabs}}) = abs.(s[n, :i])
Base.getindex(s::AbstractState, n, ::Type{Val{:δ}}) = angle.(s[n, :i])
Base.getindex(s::AbstractState, n, ::Type{Val{:s}}) = s[n, :u] .* conj.(s[n, :i])
Base.getindex(s::AbstractState, n, ::Type{Val{:p}}) = real.(s[n, :s])
Base.getindex(s::AbstractState, n, ::Type{Val{:q}}) = imag.(s[n, :s])

internalindex(s, ::Colon, i) = internalindex(s, eachindex(Nodes(s), i))
internalindex(s, n::AbstractArray, i) = internalindex.(Ref(s), n, Ref(i))
internalindex(s, n::AbstractArray, i::AbstractArray) = internalindex.(Ref(s), n, i)
internalindex(s, n::Integer, sym::Symbol) = begin
    symbolindex = findfirst(sym .== internalsymbolsof(Nodes(s)[n]))
    if symbolindex === nothing
        throw(StateError("Cannot find symbol $(QuoteNode(sym)) for node $n in the state."))
    end
    internalindex(s, n, symbolindex)
end
internalindex(s, n::Integer, i::Integer) = (NetworkRHS(s).intrange.start - 1) + NetworkRHS(s).nodalintranges[n][i]

Base.getindex(s::AbstractState, n, ::Type{Val{:int}}, i) = begin
    BaseState(s)[internalindex(s, n, i)]
end
Base.getindex(s::AbstractState, n, ::Type{Val{sym}}) where sym = getindex(s, n, Val{:int}, sym)

Base.setindex!(s::AbstractState, v, n::Colon, sym::Symbol, args...) = setindex!(s, v, eachindex(Nodes(s)), sym, args...)
Base.setindex!(s::AbstractState, v, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(Nodes(s))  )
        throw(BoundsError(s, n))
    end
    setindex!(s, v,  n, Val{sym}, args...)
end
Base.setindex!(s::AbstractState, v, n, ::Type{Val{:u}}) = begin
    setindex!(BaseState(s), real(v) ,2 .* n .- 1)
    setindex!(BaseState(s), imag(v), 2 .* n)
end
Base.setindex!(s::AbstractState, v, n, ::Type{Val{:int}}, i) = begin
    BaseState(s)[internalindex(s, n, i)] = v
end
Base.setindex!(s::AbstractState, v, n, ::Type{Val{sym}}) where sym = setindex!(s, v, n, Val{:int}, sym)

convert(::Type{BaseState}, s::State) = BaseState(s.base)
convert(::Type{A}, s::State{G, V, T}) where {V, A <:AbstractVector{V}, T, G} = @>> s BaseState convert(A)
