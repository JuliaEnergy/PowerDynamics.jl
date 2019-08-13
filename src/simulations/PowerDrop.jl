using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, reinit!, savevalues!
using Setfield

@Base.kwdef struct PowerDrop
    node_number
    fraction
    t_fault
    t_clearing
end

function (pd::PowerDrop)(powergrid)
    #TODO: check if node type supported for power drop
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    node_for_drop = @set node_for_drop.P *= pd.fraction
    node_list_power_drop[pd.node_number] = node_for_drop
    PowerGrid(powergrid.graph, node_list_power_drop, powergrid.lines)
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


export PowerDrop
export simulate
