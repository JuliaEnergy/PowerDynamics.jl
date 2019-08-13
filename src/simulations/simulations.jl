using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, reinit!, savevalues!
using LightGraphs: rem_edge!, edges, Edge
using Setfield
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
    fraction
    t_fault
    t_clearing
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
    #TODO: check if node type supported for power drop
    node_list_power_drop = copy(powergrid.nodes)
    @set node_list_power_drop[pd.node_number].P = node_list_power_drop[pd.node_number].P * pd.fraction
    PowerGrid(powergrid.graph, node_list_power_drop, powergrid.lines)
end

const iipfunc = true # is in-place function

function simulate(p::Perturbation, powergrid, x0; timespan)
    solve(powergrid, p(x0), timespan);
end

function simulate(lf::LineFault, powergrid, x0; timespan)
    solve(lf(powergrid), x0, timespan);
end

function simulate(pd::PowerDrop, powergrid, x0, timespan)
    @assert first(timespan) <= pd.t_fault "fault cannot begin in the past"
    @assert pd.t_clearing <= last(timespan) "fault cannot end in the future"

    problem = ODEProblem{iipfunc}(rhs(pd(powergrid)), x0.vec, (first(timespan), pd.t_clearing))
    fault_integrator = init(problem, Rodas4(autodiff=false))

    reinit!(fault_integrator, x0.vec, t0=pd.t_fault, tf=pd.t_clearing, erase_sol=false)
    savevalues!(fault_integrator)
    solve!(fault_integrator)

    problem = ODEProblem{iipfunc}(rhs(powergrid), fault_integrator.u, (pd.t_clearing, last(timespan)))
    integrator = init(problem, Rodas4(autodiff=false))

    # Now the trick: copy solution object to new integrator and
    # make sure the counters are updated, otherwise sol is overwritten in the
    # next step.
    integrator.sol = fault_integrator.sol
    integrator.saveiter = fault_integrator.saveiter
    integrator.saveiter_dense = fault_integrator.saveiter_dense
    integrator.success_iter = fault_integrator.success_iter

    solve!(integrator)

    return PowerGridSolution(integrator.sol, powergrid)


end

simulate(pd::PowerDrop, powergrid, x0) = simulate(pd, powergrid, x0, (0., pd.t_clearing))


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
export SinglePhaseShortCircuitToGround
export simulate
