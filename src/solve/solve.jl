using DifferentialEquations
using DiffEqBase: DEProblem
import DiffEqBase: solve

"""
    function solve(pg::PowerGrid, x0, timespan)

Solve a power grid `pg` (of type [`PowerDynBase.PowerGrid`](@ref)) starting at `x0` for a `timespan`, using DifferentialEquations.jl in the back.
The correct solvers are automatically chosen.
"""
const iipfunc = true # is in-place function

function solve(pg::PowerGrid, x0, timespan)
    problem = ODEProblem{iipfunc}(pg.network_dynamics,x0,timespan)
    solution = solve(problem, Rodas4(autodiff=false), force_dtmin=true)
    PowerGridSolution(solution, pg)
end
