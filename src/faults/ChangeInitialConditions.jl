using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

@doc """
```Julia
ChangeInitialConditions(;node,var,f)
```
# Keyword Arguments
- `node`: number  of the node
- `var`: symbol of the variable to be perturbated
- `f`: function for mapping the variable x to the perturbated value"""
struct ChangeInitialConditions
    node
    var
    f
end

ChangeInitialConditions(;node,var,f)=ChangeInitialConditions(node,var,f)

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

function (p::ChangeInitialConditions)(op)
    x0 = copy(op)
    old = x0[p.node, p.var]
    new = p.f(old)
    x0[p.node, p.var] = new
    x0
end


const iipfunc = true # is in-place function

"""
```Julia
simulate(p::ChangeInitialConditions, powergrid, x0; timespan)
```
Simulates a [`ChangeInitialConditions`](@ref)
"""
function simulate(p::ChangeInitialConditions, powergrid, x0::State, timespan; solve_kwargs...)
    solve(powergrid, p(x0), timespan; solve_kwargs...);
end


function solve(pg::PowerGrid, x0, timespan; solve_kwargs...)
    problem = ODEProblem{iipfunc}(rhs(pg),x0.vec,timespan)
    solution = solve(problem, Rodas4(); solve_kwargs...)
    PowerGridSolution(solution, pg)
end

export Inc
export Dec
export ChangeInitialConditions
export simulate