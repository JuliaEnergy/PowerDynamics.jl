using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!
using Setfield

@Base.kwdef struct PowerDrop
    node_number
    fraction
    tspan_fault
end

function (pd::PowerDrop)(powergrid)
    #TODO: check if node type supported for power drop
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    node_for_drop = @set node_for_drop.P *= pd.fraction
    node_list_power_drop[pd.node_number] = node_for_drop
    PowerGrid(node_list_power_drop, powergrid.lines)
end

function simulate(pd::PowerDrop, powergrid, x0, timespan)
    @assert first(timespan) <= pd.tspan_fault[1] "fault cannot begin in the past"
    @assert pd.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    problem = ODEProblem{true}(rhs(powergrid), x0.vec, timespan)
    integrator = init(problem, Rodas4(autodiff=false))

    step!(integrator, pd.tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs(pd(powergrid))

    step!(integrator, pd.tspan_fault[2], true)

    # update integrator, clear error
    integrator.f = rhs(powergrid)

    solve!(integrator)

    return PowerGridSolution(integrator.sol, powergrid)
end

function simulate2(pd::PowerDrop, powergrid, x0, timespan)
    @assert first(timespan) <= pd.tspan_fault[1] "fault cannot begin in the past"
    @assert pd.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    sol1 = solve(powergrid, x0, (timespan[1], pd.tspan_fault[1]))

    final_state1 = sol1(:final)
    powergrid_drop = pd(powergrid)
    sol2 = solve(powergrid_drop, final_state1, (pd.tspan_fault[1], pd.tspan_fault[2]))

    final_state2 = sol2(:final)
    # solve after clearance of the fault
    sol3 = solve(powergrid, final_state2, (pd.tspan_fault[2], last(timespan)))

    return (sol1, sol2, sol3)
end

export PowerDrop
export simulate, simulate2
