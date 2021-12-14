
function filter_lines(lines::OrderedDict, line_name)
    OrderedDict(k => v for (k, v) in lines if k != line_name)
end

function filter_lines(lines::Array, line_name)
    copy(lines[1:end.!=line_name])
end


@doc """
```Julia
LineFailure(;line_name,tspan_fault)
```
The arguments `line_name` and `tspan_fault` specify the line that should be disconnected from the grid and the respective fault duration.
For a list of lines the line_name is the index and for an OrderedDict it is the key of the line.
"""
struct LineFailure <: AbstractPerturbation
    line_name::Any
    tspan_fault::Any
    LineFailure(; line_name = line_name, tspan_fault = tspan_fault) =
        new(line_name, tspan_fault)
end

function (lf::LineFailure)(powergrid)
    filtered_lines = filter_lines(powergrid.lines, lf.line_name)
    PowerGrid(powergrid.nodes, filtered_lines)
end


export LineFailure
export simulate
