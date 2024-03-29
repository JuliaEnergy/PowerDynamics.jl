using DiffEqBase: AbstractTimeseriesSolution
using SciMLBase
using Lazy: @>>
using RecipesBase: @recipe, RecipesBase
using NetworkDynamics


"""
    struct PowerGridSolution
        dqsol::AbstractTimeseriesSolution
        powergrid::PowerGrid
    end
The data structure interfacing to the solution of the differential equations of a power grid.
Normally, it is not created by hand but return from `PowerDynSolve.solve`.
# Accessing the solution in a similar interface as [`State`](@ref).
For some grid solution `sol`, one can access the variables as
```julia
sol(t, n, s)
```
where `t` is the time (either float or array),
`n` the node number(s) (either integer, array, range (e.g. 2:3) or colon (:, for all nodes)),
and `s` is the symbol represnting the chosen value.
`s` can be either: `:v`, `:φ`, `:i`, `:iabs`, `:δ`, `:s`, `:p`, `:q`, or the symbol of the internal variable of the nodes.
The meaning of the symbols derives from the conventions of PowerDynamics.jl.
Finally, one can access the `a`-th internal variable of a node by using `sol(t, n, :int, a)`.
# Interfacing the Plots.jl library via plotting recipes, that follow similar instructions as the direct access to the solution.
For some grid solution `sol`, one plot variables of the solution asin
```julia
using Plots
plot(sol, n, s, plots_kwargs...)
```
where `n` and `s` are as in the accessing of plots, and `plots_kwargs` are the keyword arguments for Plots.jl.
"""
struct PowerGridSolution
    dqsol::AbstractTimeseriesSolution
    powergrid::PowerGrid
    function PowerGridSolution(dqsol::AbstractTimeseriesSolution, powergrid::PowerGrid)
        if !SciMLBase.successful_retcode(dqsol.retcode)
            throw(GridSolutionError("unsuccesful, return code is $(dqsol.retcode)"))
        end
        new(dqsol, powergrid)
    end
end
TimeSeries(sol::PowerGridSolution) = sol.dqsol
tspan(sol::PowerGridSolution) = (TimeSeries(sol).t[1], TimeSeries(sol).t[end])
tspan(sol::PowerGridSolution, tres) = range(TimeSeries(sol).t[1], stop=TimeSeries(sol).t[end], length=tres)

(sol::PowerGridSolution)(sym::Symbol) = sol(Val{sym})
(sol::PowerGridSolution)(::Type{Val{:initial}}) = sol(tspan(sol)[1])
(sol::PowerGridSolution)(::Type{Val{:final}}) = sol(tspan(sol)[2])


(sol::PowerGridSolution)(t::Number) = begin
    State(sol.powergrid, convert(Array, TimeSeries(sol)(t)))
end

(sol::PowerGridSolution)(t, ::Colon, sym::Symbol, args...; kwargs...) = sol(t, collect(keys(sol.powergrid.nodes)), sym, args...; kwargs...)
(sol::PowerGridSolution)(t, n, sym::Symbol, args...) = begin
    if ~all( 1 .<= n .<= length(sol.powergrid.nodes) )
        throw(StateError("Index: $n is outside of range of nodes."))
    else
        sol(t, n, Val{sym}, args...)
    end
end

(sol::PowerGridSolution)(t, n::String, sym::Symbol, args...) = begin
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=findfirst(x->x==n, bus_array)
    if (ni === nothing)
        throw(StateError("Key: $n is not in bus dictionary."))
    else
        sol(t, n, Val{sym}, args...)
    end
end

(sol::PowerGridSolution)(t, n::Array, sym::Symbol, args...) = begin
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=[findfirst(x->x==nx, bus_array) for nx in n]
    if any(ni.===nothing)||(~all( 1 .<= ni .<= length(sol.powergrid.nodes) ))
        throw(StateError("Array: $n is not in bus dictionary."))
    else
        sol(t, n, Val{sym}, args...)
    end
end

(sol::PowerGridSolution)(t::Number, n::Number, ::Type{Val{:u}}) = begin
    u_real = TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_r))
    u_imag = TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_i))
    u_real + im * u_imag
end

(sol::PowerGridSolution)(t::Number, n::String, ::Type{Val{:u}}) = begin
    u_real = TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_r))
    u_imag = TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_i))
    u_real + im * u_imag
end
(sol::PowerGridSolution)(t, n::AbstractArray, ::Type{Val{:u}}) = begin
    u_real = @>> TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_r)) convert(Array)
    u_imag = @>> TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_i)) convert(Array)
    u_real .+ im .* u_imag
end
(sol::PowerGridSolution)(t, n::String, ::Type{Val{:u}}) = begin
    u_real = @>> TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_r)) convert(Array)
    u_imag = @>> TimeSeries(sol)(t, idxs= variable_index(sol.powergrid.nodes, n, :u_i)) convert(Array)
    u_real .+ im .* u_imag
end
(sol::PowerGridSolution)(t, n, ::Type{Val{:v}}) = sol(t, n, :u) .|> abs
(sol::PowerGridSolution)(t, n, ::Type{Val{:φ}}) = sol(t, n, :u) .|> angle

(sol::PowerGridSolution)(t, n, ::Type{Val{:i}}) = get_current(sol, t, n)
(sol::PowerGridSolution)(t::Number, n, ::Type{Val{:i}}) = get_current(sol, t, n)

(sol::PowerGridSolution)(t, n, ::Type{Val{:iabs}}) = sol(t, n, :i) .|> abs
(sol::PowerGridSolution)(t, n, ::Type{Val{:δ}}) = sol(t, n, :i) .|> angle
(sol::PowerGridSolution)(t, n, ::Type{Val{:s}}) = sol(t, n, :u) .* conj.(sol(t, n, :i))
(sol::PowerGridSolution)(t, n, ::Type{Val{:p}}) = sol(t, n, :s) .|> real
(sol::PowerGridSolution)(t, n, ::Type{Val{:q}}) = sol(t, n, :s) .|> imag
(sol::PowerGridSolution)(t, n, ::Type{Val{:int}}, i) = @>> TimeSeries(sol)(t, idxs=variable_index(sol.powergrid.nodes, n, i)) convert(Array)
(sol::PowerGridSolution)(t::Number, n::String, ::Type{Val{:int}}, i::S) where {S <: Union{Number,Symbol}} = TimeSeries(sol)(t, idxs=variable_index(sol.powergrid.nodes, n, i))
(sol::PowerGridSolution)(t::Number, n::Number, ::Type{Val{:int}}, i::S) where {S <: Union{Number,Symbol}} = TimeSeries(sol)(t, idxs=variable_index(sol.powergrid.nodes, n, i))
(sol::PowerGridSolution)(t, n, ::Type{Val{sym}}) where sym = sol(t, n, Val{:int}, sym)

variable_index(nodes, n::AbstractArray, s) = map(n -> variable_index(nodes, n, s), n)

startindex(nodes, n::AbstractArray) = map(n -> startindex(nodes, n), n)

# current for a timeseries t
get_current(sol, t, n) = begin
    vertices = map(construct_vertex, collect(values(sol.powergrid.nodes)))
    edges = map(construct_edge, collect(values(sol.powergrid.lines)))
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=findfirst(x->x==n, bus_array)
    nd = network_dynamics(vertices, edges, sol.powergrid.graph)
    xt = sol.dqsol(t)
    hcat([get_current_internal(nd, x, ni) for x in xt]...)
end

# current for a single point in time
get_current(sol, t::Number, n) = begin
    vertices = map(construct_vertex, collect(values(sol.powergrid.nodes)))
    edges = map(construct_edge, collect(values(sol.powergrid.lines)))
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=findfirst(x->x==n, bus_array)
    nd = network_dynamics(vertices,edges,sol.powergrid.graph)
    x = sol.dqsol(t)
    gd = nd(x, nothing, 0.0, GetGD)
    total_current(get_dst_edges(gd, ni))
end

# current for a single point in time and for array of nodes
get_current(sol, t::Number, n::Array) = begin
    vertices = map(construct_vertex, collect(values(sol.powergrid.nodes)))
    edges = map(construct_edge, collect(values(sol.powergrid.lines)))
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=[findfirst(x->x==nx, bus_array) for nx in n]
    nd = network_dynamics(vertices,edges,sol.powergrid.graph)
    x = sol.dqsol(t)
    gd = nd(x, nothing, 0.0, GetGD)
    [total_current(get_dst_edges(gd, nx)) for nx in ni]
end

# current for array of nodes
get_current(sol, t, n::Array) = begin
    vertices = map(construct_vertex, collect(values(sol.powergrid.nodes)))
    edges = map(construct_edge, collect(values(sol.powergrid.lines)))
    bus_array=collect(keys(sol.powergrid.nodes))
    ni=[findfirst(x->x==nx, bus_array) for nx in n]
    nd = network_dynamics(vertices,edges,sol.powergrid.graph)
    xt = sol.dqsol(t)
    hcat([get_current_internal(nd, x, ni) for x in xt]...)
end


get_current_internal(nd, x, nodes) = begin
    gd = nd(x, nothing, 0.0, GetGD)
    [total_current(get_dst_edges(gd, n)) for n in collect(values(nodes))]
end

get_current_internal(nd, x, n::Number) = begin
    gd = nd(x, nothing, 0.0, GetGD)
    total_current(get_dst_edges(gd, n))
end

# define the plotting recipes

"Create the standard variable labels for power grid plots."
tslabel(sym, node) = "$(sym)$(node)"
tslabel(sym, n::AbstractArray) = tslabel.(Ref(sym), n)
tslabel(sym, node, i) = "$(sym)$(node)_$(i)"
tslabel(sym, n::AbstractArray, i) = tslabel.(Ref(sym), n, Ref(i))
"Transform the array output from DifferentialEquations.jl correctly to be used in Plots.jl's recipes."
tstransform(arr::AbstractArray{T, 1}) where T = arr
tstransform(arr::AbstractArray{T, 2}) where T = arr'

const PLOT_TTIME_RESOLUTION = 10_000 # TODO:@sabine is this high resolution really needed?

@recipe function f(sol::PowerGridSolution, ::Colon, sym::Symbol, args...)
    sol, collect(keys(sol.powergrid.nodes)), sym, args ...
end
@recipe function f(sol::PowerGridSolution, n, sym::Symbol, args...; tres = PLOT_TTIME_RESOLUTION)
    if sym == :int
        label --> tslabel(sym, n, args[1])
    else
        label --> tslabel(sym, n)
    end
    xguide --> "t"
    t = tspan(sol, tres)
    println(t)
    t, tstransform(sol(t, n, sym, args...))
end
