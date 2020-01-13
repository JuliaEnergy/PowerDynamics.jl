using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

"""
```Julia
Perturbation(;node,var,f)
```
# Keyword Arguments
- `node`: number  of the node
- `var`: symbol of the variable to be perturbated
- `f`: function for mapping the variable x to the perturbated value
"""
struct Perturbation
    node
    var
    f
end

Perturbation(;node,var,f)=Perturbation(node,var,f)

"""
```Julia
LineFault(;from,to)
```
The arguments `from` and `to` specify the line that should be disconnected from the grid.
"""
struct LineFault
    from
    to
end

LineFault(;from,to) = LineFault(from,to)

struct Inc
    val
end

function (i::Inc)(x)
    x += i.val
end

struct Dec
    val
end

function (i::Dec)(x)
    x -= i.val
end

function (lf::LineFault)(powergrid)
    @assert lf.from < lf.to "order important to remove the line from powergrid"
    filtered_lines = filter(l -> (l.from !=lf.from || l.to !=lf.to), copy(powergrid.lines))
    PowerGrid(powergrid.nodes, filtered_lines)
end

function (p::Perturbation)(op)
    x0 = copy(op)
    old = x0[p.node, p.var]
    new = p.f(old)
    x0[p.node, p.var] = new
    x0
end


const iipfunc = true # is in-place function

"""
```Julia
simulate(p::Perturbation, powergrid, x0; timespan)
```
Simulates a [`Perturbation`](@ref)
"""
function simulate(p::Perturbation, powergrid, x0; timespan)
    solve(powergrid, p(x0), timespan);
end

"""
```Julia
simulate(lf::LineFault, powergrid, x0; timespan)
```
Simulates a [`LineFault`](@ref)
"""
function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end

function solve(pg::PowerGrid, x0, timespan)
    problem = ODEProblem{iipfunc}(rhs(pg),x0,timespan)
    solution = solve(problem, Rodas4(autodiff=false))
    PowerGridSolution(solution, pg)
end

export Inc
export Dec
export LineFault
export Perturbation
export simulate
