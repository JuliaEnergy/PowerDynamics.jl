using OrdinaryDiffEq: ODEProblem, Rodas4
using LightGraphs: rem_edge!
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

@Base.kwdef struct PowerDrop
    percent
    node_number
    t1
    t2
    t3
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

function (pd::PowerDrop)(powergrid)
    node_list_power_drop = copy(powergrid.nodes)
    n = powergrid.nodes[pd.node_number]
    node_list_power_drop[pd.node_number] = SwingEqLVS(
        H=n.H,
        P=n.P*pd.percent,
        D=n.D,
        Ω=n.Ω,
        Γ=n.Γ,
        V=n.V
    )
    PowerGrid(powergrid.graph, node_list_power_drop, powergrid.lines)
end

function simulate(p::Perturbation, powergrid, x0; timespan)
    solve(powergrid, p(x0), timespan);
end

function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end

function simulate(pd::PowerDrop, powergrid, x0)
    sol1 = solve(powergrid, x0, pd.t1)
    final_state1 = sol1(:final)

    g_power_reduction = pd(powergrid)
    sol2 = solve(g_power_reduction, State(g_power_reduction, final_state1.vec), pd.t2)
    final_state2 = sol2(:final)

    sol3 = solve(powergrid, State(powergrid, final_state2.vec), pd.t3)
    sol3
end

const iipfunc = true # is in-place function

function solve(pg::PowerGrid, x0, timespan)
    problem = ODEProblem{iipfunc}(rhs(pg),x0.vec,timespan)
    solution = solve(problem, Rodas4(autodiff=false))
    PowerGridSolution(solution, pg)
end

export Inc
export Dec
export LineFault
export Perturbation
export PowerDrop
export simulate
