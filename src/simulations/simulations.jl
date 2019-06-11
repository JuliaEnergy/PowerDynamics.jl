using DifferentialEquations
using DiffEqBase: DEProblem
import DiffEqBase: solve

@Base.kwdef struct Perturbation
    node
    var
    f
end

@Base.kwdef struct LineFault
    from
    to
end

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
    graph_copy = copy(powergrid.graph)
    rem_edge!(graph_copy, lf.from, lf.to)
    PowerGrid(graph_copy)
end

function (p::Perturbation)(op)
    x0 = copy(op)
    old = x0[p.node, p.var]
    new = p.f(old)
    x0[p.node, p.var] = new
    x0
end

function simulate(p::Perturbation, powergrid, x0; timespan)
    solve(powergrid, p(x0), timespan);
end

function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end

const iipfunc = true # is in-place function

function solve(pg::PowerGrid, x0, timespan)
    problem = ODEProblem{iipfunc}(ode_function(pg),x0.vec,timespan)
    solution = solve(problem, Rodas4(autodiff=false))
    PowerGridSolution(solution, pg)
end

export Inc
export Dec
export LineFault
export Perturbation
export simulate
