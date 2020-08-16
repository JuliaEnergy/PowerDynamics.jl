using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

"""
```Julia
NodePerturbation(;node,var,f,timespan)
```
# Keyword Arguments
- `node`: number  of the node
- `var`: symbol of the variable to be perturbated
- `f`: function for mapping the variable x to the perturbated value
- `timespan` time span for which perturbation is present
"""
abstract type NodePerturbation end

"Error to be thrown if something goes wrong during power perturbation"
struct NodePerturbationError <: PowerDynamicsError
    msg::String
end

function mapField(powergrid, np, f)
    nodes = copy(powergrid.nodes)
    np_node = powergrid.nodes[np.node]

    # add shunt to the node component at the fault location
    if !(hasproperty(np_node, np.var))
        throw(NodePerturbationError("Node number: $(np.node) must have a variable called $(np.var)."))
    end
    lens = Setfield.compose(Setfield.PropertyLens{np.var}())
    node_for_perturbation = Setfield.modify(f, np_node, lens)
    nodes[np.node] = node_for_perturbation
    PowerGrid(nodes, powergrid.lines)
end

"""
```Julia
simulate(no::NodePerturbation, powergrid, x1, timespan)
```
Simulates a [`NodePerturbation`](@ref)
"""
function simulate(np::NodePerturbation, powergrid, x1, timespan)
    @assert first(timespan) <= np.tspan_fault[1] "fault cannot begin in the past"
    @assert np.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    np_powergrid = np(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)

    function errorState(integrator)
        sol1 = integrator.sol
        x2 = find_valid_initial_condition(np_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(np_powergrid)
        integrator.u = x2
    end

    function regularState(integrator)
        sol2 = integrator.sol
        x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(powergrid)
        integrator.u = x3
    end

    t1 = np.tspan_fault[1]
    t2 = np.tspan_fault[2]

    cb1 = DiscreteCallback(((u,t,integrator) -> t in np.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in np.tspan_fault[2]), regularState)
    sol = solve(problem, Rodas4(), force_dtmin = true, callback = CallbackSet(cb1, cb2), tstops=[t1, t2])
    return PowerGridSolution(sol, powergrid)
end

simulate(np::NodePerturbation, op::State, timespan) = simulate(np, op.grid, op.vec, timespan)

export NodePerturbation
export NodePerturbationError
export simulate
export mapField
