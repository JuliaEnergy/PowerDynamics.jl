using DifferentialEquations: ODEProblem, Rodas4
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
    node_number
    percent
    t_prefault
    t_fault
    t_postfault
end

@Base.kwdef struct SinglePhaseShortCircuitToGround
    #ToDo use from/to and automatically find the line
    # what is the indexing of lines?
    line_number
    line_fraction
    short_circuit_admittance
    t_prefault
    t_fault
    t_postfault
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

function (pd::PowerDrop)(powergrid)
    # ToDo type SlackAlgebraic has no field P
    # PROBLEM: node types are immutable. how to create a copy?
    # fieldnames(::NodeType) is not defined
    node_list_power_drop = copy(powergrid.nodes)
    node_list_power_drop[pd.node_number].P *= pd.percent
    PowerGrid(powergrid.graph, node_list_power_drop, powergrid.lines)
end

function (sc::SinglePhaseShortCircuitToGround)(powergrid)
    line_list_power_drop = copy(powergrid.lines)
    healthy_line = line_list_power_drop[sc.line_number]

    X = inv(healthy_line.Y)
    Xg = inv(sc.short_circuit_admittance)

    Xprime = X + (1. - sc.line_fraction) * sc.line_fraction * X * X / Xg
    Xl_shunt = inv(1. - sc.line_fraction) * Xg + sc.line_fraction * X
    Xr_shunt = inv(sc.line_fraction) * Xg + (1. - sc.line_fraction) * X

    line_list_power_drop[sc.line_number] = PiModelLine(y=inv(Xprime), y_shunt_km=inv(Xl_shunt), y_shunt_mk=inv(Xr_shunt))

    PowerGrid(powergrid.graph, powergrid.nodes, line_list_power_drop)
end

function simulate(p::Perturbation, powergrid, x0; timespan)
    solve(powergrid, p(x0), timespan);
end

function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end

function simulate(pd::PowerDrop, powergrid, x0)
    sol1 = solve(powergrid, x0, pd.t_prefault)
    final_state1 = sol1(:final)

    g_power_reduction = pd(powergrid)
    sol2 = solve(g_power_reduction, State(g_power_reduction, final_state1.vec), pd.t_fault)
    final_state2 = sol2(:final)

    sol3 = solve(powergrid, State(powergrid, final_state2.vec), pd.t_postfault)

    CompositePowerGridSolution([sol1, sol2, sol3], [powergrid, g_power_reduction, powergrid])
end

function simulate(sc::SinglePhaseShortCircuitToGround, powergrid, x0)
    sol1 = solve(powergrid, x0, sc.t_prefault)
    final_state1 = sol1(:final)

    g_faulted_line = sc(powergrid)
    sol2 = solve(g_faulted_line, State(g_faulted_line, final_state1.vec), sc.t_fault)
    final_state2 = sol2(:final)

    sol3 = solve(powergrid, State(powergrid, final_state2.vec), sc.t_postfault)

    CompositePowerGridSolution([sol1, sol2, sol3], [powergrid, g_faulted_line, powergrid])
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
export PowerDrop
export SinglePhaseShortCircuitToGround
export simulate
