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
- `power_symbol`: parameter symbol on the node that represents Power, default is :P
"""
Base.@kwdef struct PowerPerturbation
    node_number
    fraction
    tspan_fault
    power_symbol = :P
end

"Error to be thrown if something goes wrong during power perturbation"
struct PowerPerturbationError <: PowerDynamicsError
    msg::String
end

function (pd::PowerPerturbation)(powergrid)
    mapPowerField(powergrid, pd, p -> p * pd.fraction)
end

function mapPowerField(powergrid, pd, f)
    node_list_power_drop = copy(powergrid.nodes)
    node_for_drop = node_list_power_drop[pd.node_number]
    if !(hasproperty(node_for_drop,pd.power_symbol))
        throw(PowerPerturbationError("Node number: $(pd.node_number) must have a power parameter $(pd.power_symbol)"))
    end
    lens = Setfield.PropertyLens{pd.power_symbol}()
    node_for_drop = Setfield.modify(f, node_for_drop, lens)
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

    typeStablePowerGrid = mapPowerField(powergrid, pd, p -> convert(Float64, p))
    normal_rhs = rhs(typeStablePowerGrid)

    problem = ODEProblem{true}(normal_rhs, x0.vec, timespan)
    # the fault times need to be added as t_stops values to the integrator so
    # it makes extra small steps around the discontinuities
    dt_fault = 1e-7
    t1=pd.tspan_fault[1]-dt_fault
    t2=pd.tspan_fault[1]+dt_fault
    t3=pd.tspan_fault[2]-dt_fault
    t4=pd.tspan_fault[2]+dt_fault

    integrator = init(problem, Rodas4(),tstops=[t1,t2,t3,t4])

    step!(integrator, pd.tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs(pd(typeStablePowerGrid))
    u_modified!(integrator,true)

    step!(integrator, pd.tspan_fault[2]-pd.tspan_fault[1], true)

    # update integrator, clear error
    integrator.f = normal_rhs
    u_modified!(integrator,true)

    step!(integrator, timespan[2]-pd.tspan_fault[2], true)

    solve!(integrator)

    return PowerGridSolution(integrator.sol, powergrid)
end

export PowerPerturbation
export PowerPerturbationError
export simulate
