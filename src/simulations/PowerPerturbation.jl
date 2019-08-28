using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!
using Setfield

@Base.kwdef struct PowerPerturbation
    node_number
    fraction
    tspan_fault
end

function (pd::PowerPerturbation)(powergrid)
    #TODO: check if node type supported for power drop
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    node_for_drop = @set node_for_drop.P *= pd.fraction
    node_list_power_drop[pd.node_number] = node_for_drop
    PowerGrid(node_list_power_drop, powergrid.lines)
end

# this works, but the code looks hacky. Can't we do any better?
function simulate(pd::PowerPerturbation, powergrid, x0, timespan)
    @assert first(timespan) <= pd.tspan_fault[1] "fault cannot begin in the past"
    @assert pd.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    problem = ODEProblem{true}(rhs(pd(powergrid)), x0.vec, (first(timespan), pd.tspan_fault[2]))
    fault_integrator = init(problem, Rodas4(autodiff=false))

    reinit!(fault_integrator, x0.vec, t0=pd.tspan_fault[1], tf=pd.tspan_fault[2], erase_sol=false)
    savevalues!(fault_integrator)
    solve!(fault_integrator)

    problem = ODEProblem{true}(rhs(powergrid), fault_integrator.u, (pd.tspan_fault[2], last(timespan)))
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

export PowerPerturbation
export simulate, simulate2, simulate3
