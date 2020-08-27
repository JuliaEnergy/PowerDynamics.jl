using OrdinaryDiffEq: ODEProblem, Rodas4, DiscreteCallback, CallbackSet
import DiffEqBase: solve
using Setfield

abstract type AbstractPerturbation end

function typestable_node_field_update(powergrid::PowerGrid, node, sym::Symbol, val)
    old_node = powergrid.nodes[node]

    if !(hasproperty(old_node, sym))
        throw(FieldUpdateError("Node $(node) must have a variable called $(sym)."))
    end
    
    nodes = copy(powergrid.nodes)
    sym_type = getfield(old_node, sym) |> eltype
    lens = Setfield.compose(Setfield.PropertyLens{sym}())

    new_val = try
        convert(sym_type, val)
    catch e
        throw(FieldUpdateError("The parameter $(sym) at node $(node) should be of type $(sym_type). \n $(e)"))
    end

    new_node = Setfield.set(old_node, lens, new_val)
    nodes[node] = new_node
    PowerGrid(nodes, powergrid.lines)
end

"""
```Julia
simulate(no::AbstractPerturbation, powergrid, x1, timespan)
```
Simulates a [`AbstractPerturbation`](@ref)
"""
function simulate(np::AbstractPerturbation, powergrid::PowerGrid, x1, timespan; solve_kwargs...)
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

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1, cb2), tstops=[t1, t2], solve_kwargs...)
    
    return PowerGridSolution(sol, powergrid)
end

simulate(np::AbstractPerturbation, op::State, timespan) = simulate(np, op.grid, op.vec, timespan)

"Error to be thrown if something goes wrong during power perturbation"
struct FieldUpdateError <: PowerDynamicsError
    msg::String
end

export AbstractPerturbation
export simulate
export typestable_node_field_update
export FieldUpdateError
