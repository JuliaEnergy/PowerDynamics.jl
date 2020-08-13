using OrdinaryDiffEq: ODEProblem, Rodas4, solve!, DiscreteCallback, CallbackSet
import DiffEqBase: solve
using Setfield

"""
```Julia
    NodeShortCircuit(;node_number,Y,tspan_fault)
```
# Keyword Arguments
- `node_number`: number  of the node
- `Y`: admittance of the short circuit
- `tspan_fault`: short circuit timespan
"""
Base.@kwdef struct NodeShortCircuit
    node_number
    Y
    tspan_fault
    shunt_symbol = :Y_n
end

"Error to be thrown if something goes wrong during short circuit"
struct NodeShortCircuitError <: PowerDynamicsError
    msg::String
end

function (nsc::NodeShortCircuit)(powergrid)
    nodes = copy(powergrid.nodes)
    sc_node = powergrid.nodes[nsc.node_number]

    # add shunt to the node component at the fault location
    if !(hasproperty(sc_node, nsc.shunt_symbol))
        throw(NodeShortCircuitError("Node number: $(nsc.node_number) must have a shunt field called $(nsc.shunt_symbol)."))
    end
    lens = Setfield.compose(Setfield.PropertyLens{nsc.shunt_symbol}())
    faulted_component = Setfield.set(sc_node, lens, nsc.Y)

    nodes[nsc.node_number] = faulted_component
    PowerGrid(nodes, powergrid.lines)
end

"""
```Julia
simulate(nsc::NodeShortCircuit, powergrid, x1, timespan)
```
Simulates a [`NodeShortCircuit`](@ref)
"""
function simulate(nsc::NodeShortCircuit, powergrid, x1, timespan)
    @assert first(timespan) <= nsc.tspan_fault[1] "fault cannot begin in the past"
    @assert nsc.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    nsc_powergrid = nsc(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)

    function errorState(integrator)
        sol1 = integrator.sol
        x2 = find_valid_initial_condition(nsc_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(nsc_powergrid)
        integrator.u = x2
    end

    function regularState(integrator)
        sol2 = integrator.sol
        x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(powergrid)
        integrator.u = x3
    end

    t1 = nsc.tspan_fault[1]
    t2 = nsc.tspan_fault[2]

    cb1 = DiscreteCallback(((u,t,integrator) -> t in nsc.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in nsc.tspan_fault[2]), regularState)
    sol = solve(problem, Rodas4(), force_dtmin = true, callback = CallbackSet(cb1, cb2), tstops=[t1, t2])
    return PowerGridSolution(sol, powergrid)
end

simulate(nsc::NodeShortCircuit, op::State, timespan) = simulate(nsc, op.grid, op.vec, timespan)

export NodeShortCircuit
export NodeShortCircuitError
export simulate
