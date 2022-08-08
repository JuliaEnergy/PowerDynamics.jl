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
- Use a [`PowerGrid`](@ref) instance `grid` and a properly sized
    state vector `vec` to instantiate a `State`.

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

You can access (and set) the ``k``-th variable by calling

`s[j, :var, k]`.

The variables can be also directly accessed with symbols, i.e.

`s[j, :ω]`

returns the frequency ``ω`` at node ``j``.
To find out the proper symbol, the easiest way is to look into the docs
of the corresponding node type,
check the output of symbolsof or simply look at the output of `println`:

    julia> symbolsof(SwingEq(H=2, P=3, D=4, Ω=5))
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

Base.getindex(s::State, n::Colon, sym::Symbol, args...) = [getindex(s, ni, sym, args...) for ni in collect(keys(s.grid.nodes))]
Base.getindex(s::State, n::Array, sym::Symbol, args...) = [getindex(s, ni, sym, args...) for ni in n]

Base.getindex(s::State, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(s.grid.nodes)  )
        throw(BoundsError(s, n))
    end
    getindex(s, n, Val{sym}, args...)
end

Base.getindex(s::State, n::String, sym::Symbol, args...) = begin
    bus_array=collect(keys(s.grid.nodes))
    ni=findfirst(x->x==n, bus_array) #for nx in Array([n])]
    #ni = ni[1]:ni[end]
    if ni === nothing
        throw(StateError("Node $n is not in bus dictionary."))
    end
    getindex(s, n, Val{sym}, args...)
end

Base.getindex(s::State, n, ::Type{Val{:u}}) = getindex(s, n, :u_r) + im .* getindex(s, n, :u_i)

Base.getindex(s::State, n, ::Type{Val{:v}}) = abs.(s[n, :u])
Base.getindex(s::State, n, ::Type{Val{:φ}}) = angle.(s[n, :u])
Base.getindex(s::State, n, ::Type{Val{:i}}) = get_current.(Ref(s), n)
Base.getindex(s::State, n, ::Type{Val{:iabs}}) = abs.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:δ}}) = angle.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:s}}) = s[n, :u] .* conj.(s[n, :i])
Base.getindex(s::State, n, ::Type{Val{:p}}) = real.(s[n, :s])
Base.getindex(s::State, n, ::Type{Val{:q}}) = imag.(s[n, :s])

internalindex(s, ::Colon, i) = internalindex(s, eachindex(s.grid.nodes, i))
internalindex(s, n::AbstractArray, i) = internalindex.(Ref(s), n, Ref(i))
internalindex(s, n::AbstractArray, i::AbstractArray) = internalindex.(Ref(s), n, i)
internalindex(s, n::Integer, sym::Symbol) = variable_index(s.grid.nodes, n, sym)
internalindex(s, n::Integer, i) = variable_index(s.grid.nodes, n, i)
internalindex(s, n::String, sym::Symbol) = variable_index(s.grid.nodes, n, sym)
internalindex(s, n::String, i) = variable_index(s.grid.nodes, n, i)
internalindex(s, n::String, i::AbstractArray) = internalindex.(Ref(s), n, i)

variable_index(nodes, n::String, s::Symbol) = begin
    first_idx = findfirst(ns -> ns == s, symbolsof(nodes[n]))
    if first_idx === nothing
        throw(StateError("Variable: $s not defined for node: $n"))
    else
        startindex(nodes, n) + first_idx
    end
end

variable_index(nodes, n::String, i) = begin
    num_vars = dimension(nodes[n])
    if i <= num_vars
        startindex(nodes, n) + i
    else
        throw(BoundsError("Variable index: $i not supported for node: $(nodes[n])"))
    end
end

variable_index(nodes, n::Integer, s::Symbol) = begin
    first_idx = findfirst(ns -> ns == s, symbolsof(nodes[n]))
    if first_idx === nothing
        throw(StateError("Variable: $s not defined for node: $n"))
    else
        startindex(nodes, n) + first_idx
    end
end

variable_index(nodes, n, i) = begin
    num_vars = dimension(nodes[n]) #    #num_vars = sum([dimension(nodes[ni]) for ni in n])
    if i <= num_vars
        startindex(nodes, n) + i
    else
        throw(BoundsError("Variable index: $i not supported for node: $(nodes[n])"))
    end
end

@views startindex(nodes, n) = begin
    bus_array=collect(keys(nodes))
    values_array=collect(values(nodes))
    ni=findfirst(x->x==n, bus_array)
    if ni == 1
        0
    else
        sum(map(node -> dimension(node), values_array[1:ni-1]))
    end
end

Base.getindex(s::State, n, ::Type{Val{:var}}, i) = begin
    s[internalindex(s, n, i)]
end
Base.getindex(s::State, n::String, ::Type{Val{:var}}, i) = begin
    s[internalindex(s, n, i)]
end
Base.getindex(s::State, n, ::Type{Val{sym}}) where sym = getindex(s, n, Val{:var}, sym)


Base.setindex!(s::State, v, n::Colon, sym::Symbol, args...) = setindex!(s, v, collect(keys(s.grid.nodes)), sym, args...)

Base.setindex!(s::State, v, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(s.grid.nodes)  )
        throw(BoundsError(s, n))
    end
    setindex!(s, v,  n, Val{sym}, args...)
end

Base.setindex!(s::State, v, n::String, sym::Symbol, args...) = begin
    bus_array=collect(keys(s.grid.nodes))
    ni=findfirst(x->x==n, bus_array)
    if (ni === nothing)
        throw(StateError("Node with key: $n is not defined."))
    end
    setindex!(s, v,  n, Val{sym}, args...)
end

Base.setindex!(s::State, v, n::Array, sym::Symbol, args...) = begin
    bus_array=collect(keys(s.grid.nodes))
    ni=[findfirst(x->x==nx, bus_array) for nx in n]
    if ~all( 1 .<= ni .<= length(s.grid.nodes)  )
        throw(BoundsError(s, ni))
    end
    setindex!(s, v,  n, Val{sym}, args...)
end

Base.setindex!(s::State, v, n, ::Type{Val{:u}}) = begin
    setindex!(s, real(v), variable_index(s.grid.nodes, n, :u_r))
    setindex!(s, imag(v), variable_index(s.grid.nodes, n, :u_i))
end
Base.setindex!(s::State, v, n, ::Type{Val{:v}}) = begin
    u = s[n, :u]
    s[n, :u] = v.*u./abs.(u)
end
Base.setindex!(s::State, v, n, ::Type{Val{:φ}}) = begin
    u = s[n, :u]
    s[n, :u] = abs.(u).*exp.(im.*v)
end
Base.setindex!(s::State, v, n, ::Type{Val{:var}}, i) = begin
    s[internalindex(s, n, i)] = v
end
Base.setindex!(s::State, v, n, ::Type{Val{sym}}) where sym = setindex!(s, v, n, Val{:var}, sym)

function Base.:+(s1::State, s2::State)
    @assert s1.grid == s2.grid
    State(s1.grid, s1.vec .+ s2.vec)
end
function Base.:-(s1::State, s2::State)
    @assert s1.grid == s2.grid
    State(s1.grid, s1.vec .- s2.vec)
end

function Base.:(==)(s1::State, s2::State)
    all(s1.vec .== s2.vec)
end
function Base.:≈(s1::State, s2::State)
    all(s1.vec .≈ s2.vec)
end

Base.:-(s::State) = State(s.grid, .- s.vec)
Base.:*(k, s::State) = State(s.grid, k.*s.vec)
Base.:*(s::State, k) = k*s # commutivity
Base.:/(s::State, k) = (1/k)*s

get_current(state, n) = begin
    vertices = map(construct_vertex, state.grid.nodes)
    edges = map(construct_edge, state.grid.lines)
    nd = network_dynamics(vertices,edges,state.grid.graph)
    gd = nd(state.vec, nothing, 0.0, GetGD)
    total_current(get_dst_edges(gd, n))
end

get_current(state, n::String) = begin
    busnames_array = collect(keys(state.grid.nodes))
    node_array=collect(values(state.grid.nodes))
    vertices = map(construct_vertex, node_array)
    bus_id = findfirst(x->x==n, busnames_array)
    edges_array = collect(values(state.grid.lines))
    edges = map(construct_edge, edges_array)
    nd = network_dynamics(vertices,edges,state.grid.graph)
    gd = nd(state.vec, nothing, 0.0, GetGD)
    total_current(get_dst_edges(gd, bus_id))
end
