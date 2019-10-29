using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!
using Setfield

"""
```Julia
PowerPerturbation(;node_number,fraction,tspan_fault)
```
# Keyword Arguments
- `node_number`: number  of the node
- `fraction`: percentage factor to be applied to the active power P
- `tspan_fault`: PowerPerturbation timespan
"""
struct PowerPerturbation
    node_number
    fraction
    tspan_fault
end
PowerPerturbation(;node_number,fraction,tspan_fault) = PowerPerturbation(node_number,fraction,tspan_fault)

function (pd::PowerPerturbation)(powergrid)
    #TODO: check if node type supported for power drop
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    node_for_drop = @set node_for_drop.P *= pd.fraction
    node_list_power_drop[pd.node_number] = node_for_drop
    PowerGrid(node_list_power_drop, powergrid.lines)
end

"""
```Julia
simulate(nsc::PowerPerturbation, powergrid, x0, timespan)
```
Simulates a [`PowerPerturbation`](@ref)
"""
function simulate(pd::PowerPerturbation, powergrid, x0, timespan)
    @assert first(timespan) <= pd.tspan_fault[1] "fault cannot begin in the past"
    @assert pd.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    problem = ODEProblem{true}(rhs(powergrid), x0.vec, timespan)
    integrator = init(problem, Rodas4(autodiff=false))

    step!(integrator, pd.tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs(pd(powergrid))
    u_modified!(integrator,true)

    step!(integrator, pd.tspan_fault[2], true)

    # update integrator, clear error
    integrator.f = rhs(powergrid)
    u_modified!(integrator,true)

    step!(integrator, timespan[2]-pd.tspan_fault[2], true)


    solve!(integrator)

    return PowerGridSolution(integrator.sol, powergrid)
end

function simulate2(pd::PowerPerturbation, powergrid, x1, timespan)
    @assert first(timespan) <= pd.tspan_fault[1] "fault cannot begin in the past"
    @assert pd.tspan_fault[2] <= last(timespan) "fault cannot end in the future"
    pd_powergrid = pd(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)

    function errorState(integrator)
        sol1 = integrator.sol
        x2 = operationpoint = find_operationpoint(pd_powergrid)
        integrator.f = rhs(pd_powergrid)
        integrator.u = x2
    end

    function regularState(integrator)
        sol2 = integrator.sol
        x3 = find_operationpoint(pd_powergrid) # Jump the state to be valid for the new system.
        integrator.f = rhs(powergrid)
        integrator.u = x3
    end

    cb1 = DiscreteCallback(((u,t,integrator) -> t in pd.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in pd.tspan_fault[2]), regularState)
    sol = solve(problem, Rodas4(autodiff=false), callback = CallbackSet(cb1, cb2), tstops=[pd.tspan_fault[1];pd.tspan_fault[2]])
    return PowerGridSolution(sol, powergrid)
end

export PowerPerturbation
export simulate2
export simulate
