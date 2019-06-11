# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

import Base: convert


"""
```Julia

    State(base; t=nothing)
    State(grid, vec; t=nothing)

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
an array.
There are plenty of convenient ways to access its contents at a node ``j``
by using a particular symbol:

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
of the corresponding [`AbstractNodeParameters`](@ref) subtype,
check the output of [`internalsymbolsof`](@ref)
or simply look at the output of `println`:

    julia> internalsymbolsof(SwingEq(H=2, P=3, D=4, Ω=5))
    1-element Array{Symbol,1}:
     :ω

    julia> println(SwingEq(H=2, P=3, D=4, Ω=5))
    SwingEq[:ω](H=2, P=3, D=4, Ω=5)

"""
struct State{G <: PowerGrid, V}
    grid::G
    vec::AbstractVector{V}
    function State{G, V}(grid::G, vec::AbstractVector{V}) where {G <: PowerGrid, V}
        @assert systemsize(grid) == length(vec)
        new{G, V}(grid, vec)
    end
end
State(grid::G, vec::AbstractVector{V}) where {G <: PowerGrid, V} = State{G, V}(grid, vec)
State(grid::G, V::Type) where {G <: PowerGrid} = State{G, V}(grid)
State(grid::G) where {G <: PowerGrid} = State(grid, Float64)

Base.size(s::State) = size(s.vec)
Base.getindex(s::State, n) = getindex(s.vec, n)
Base.setindex!(s::State, v,  n) = setindex!(s.vec, v, n)
Base.copy(s::State) = State(s.grid, deepcopy(s.vec)) # grid should not be copied

convert(::Type{A}, s::State) where {A<:AbstractArray} = convert(A, s.vec)

Base.getindex(s::State, n::Colon, sym::Symbol, args...) = getindex(s, eachindex(s.grid.nodes), sym, args...)
Base.getindex(s::State, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(s.grid.nodes)  )
        throw(BoundsError(s, n))
    end
    getindex(s, n, Val{sym}, args...)
end
Base.getindex(s::State, n, ::Type{Val{:u}}) = getindex(s, 2 .* n .- 1) + im .* getindex(s, 2 .* n)
Base.getindex(s::State, n, ::Type{Val{:v}}) = abs.(s[n, :u])
Base.getindex(s::State, n, ::Type{Val{:φ}}) = angle.(s[n, :u])
#Base.getindex(s::State, n, ::Type{Val{:i}}) = getindex(AdmittanceLaplacian(s) * s[:, :u], n)
Base.getindex(s::State, n, ::Type{Val{:iabs}}) = abs.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:δ}}) = angle.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:s}}) = s[n, :u] .* conj.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:p}}) = real.(s[n, :s])
Base.getindex(s::State, n, ::Type{Val{:q}}) = imag.(s[n, :s])

internalindex(s, ::Colon, i) = internalindex(s, eachindex(Nodes(s), i))
internalindex(s, n::AbstractArray, i) = internalindex.(Ref(s), n, Ref(i))
internalindex(s, n::AbstractArray, i::AbstractArray) = internalindex.(Ref(s), n, i)
internalindex(s, n::Integer, sym::Symbol) = variable_index(s.grid.nodes, n, sym)

variable_index(nodes, n, s::Symbol) = startindex(nodes, n) + findfirst(ns -> ns == s, symbolsof(nodes[n]))

@views startindex(nodes, n) = begin
    if n == 1
        0
    else
        sum(map(node -> dimension(node), nodes[1:n-1]))
    end
end

Base.getindex(s::State, n, ::Type{Val{:int}}, i) = begin
    s[internalindex(s, n, i)]
end
Base.getindex(s::State, n, ::Type{Val{sym}}) where sym = getindex(s, n, Val{:int}, sym)


Base.setindex!(s::State, v, n::Colon, sym::Symbol, args...) = setindex!(s, v, eachindex(s.grid.nodes), sym, args...)
Base.setindex!(s::State, v, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(s.grid.nodes)  )
        throw(BoundsError(s, n))
    end
    setindex!(s, v,  n, Val{sym}, args...)
end
Base.setindex!(s::State, v, n, ::Type{Val{:u}}) = begin
    setindex!(s, real(v) ,2 .* n .- 1)
    setindex!(s, imag(v), 2 .* n)
end
Base.setindex!(s::State, v, n, ::Type{Val{:v}}) = begin
    u = s[n, :u]
    s[n, :u] = v.*u./abs.(u)
end
Base.setindex!(s::State, v, n, ::Type{Val{:int}}, i) = begin
    s[internalindex(s, n, i)] = v
end
Base.setindex!(s::State, v, n, ::Type{Val{sym}}) where sym = setindex!(s, v, n, Val{:int}, sym)

function Base.:(==)(s1::State, s2::State)
    all(s1.vec .== s2.vec)
end
function Base.:≈(s1::State, s2::State)
    all(s1.vec .≈ s2.vec)
end
