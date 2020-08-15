
"""
```Julia
LineFailure(;from,to)
```
The arguments `from` and `to` specify the line that should be disconnected from the grid.
"""
struct LineFailure(;from,to)
    from
    to
end

function (lf::LineFailure)(powergrid)
    @assert lf.from < lf.to "order important to remove the line from powergrid"
    filtered_lines = filter(l -> (l.from !=lf.from || l.to !=lf.to), copy(powergrid.lines))
    PowerGrid(powergrid.nodes, filtered_lines)
end

"""
```Julia
simulate(lf::LineFailure, powergrid, x0; timespan)
```
Simulates a [`LineFailure`](@ref)
"""
function simulate(lf::LineFailure, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end



LineFailure(;from,to) = LineFailure(from,to)
