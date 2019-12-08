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
    power_symbol
end

"Error to be thrown if something goes wrong during power perturbation"
struct PowerPerturbationError <: PowerDynamicsError
    msg::String
end

PowerPerturbation(;node_number,fraction,tspan_fault,power_symbol) = PowerPerturbation(node_number,fraction,tspan_fault,power_symbol)

function (pd::PowerPerturbation)(powergrid)
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    if !(pd.power_symbol in fieldnames(typeof(node_for_drop)))
        throw(PowerPerturbationError("Node number: $(pd.node_number) must have a power parameter: $(pd.power_symbol)"))
    end
    try
        new_P = typeof(node_for_drop.P)(node_for_drop.P * pd.fraction)
        node_for_drop = @set node_for_drop.P = new_P
        node_list_power_drop[pd.node_number] = node_for_drop
        PowerGrid(node_list_power_drop, powergrid.lines)
    catch e
        if isa(e, InexactError)
            throw(PowerPerturbationError("Node number: $(pd.node_number) must have a $(pd.power_symbol) parameter of type Float64 for PowerPerturbation, but was of type $(typeof(node_for_drop.P))"))
        else
            throw(x) #handle every other error
        end
    end
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

    step!(integrator, pd.tspan_fault[2]-pd.tspan_fault[1], true)

    # update integrator, clear error
    integrator.f = rhs(powergrid)
    u_modified!(integrator,true)

    step!(integrator, timespan[2]-pd.tspan_fault[2], true)


    solve!(integrator)

    return PowerGridSolution(integrator.sol, powergrid)
end

export PowerPerturbation
export PowerPerturbationError
export simulate
