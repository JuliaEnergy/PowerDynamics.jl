<<<<<<< HEAD
export NodeShortCircuit, simulate
using OrdinaryDiffEq: ODEProblem, Rodas4, set_u!,init,reinit!, solve!, step!, reinit!, savevalues!, u_modified!, DiscreteCallback, CallbackSet
=======
export NodeShortCircuit
using OrdinaryDiffEq: ODEProblem, Rodas4, solve!, DiscreteCallback, CallbackSet
>>>>>>> cleanup
import DiffEqBase: solve

"""
```Julia
NodeShortCircuit(;node,var,f)
```
# Keyword Arguments
- `node`: number  of the node
- `R`: resistance of the shortcircuit
- `sc_timespan`: shortcircuit timespan
"""
struct NodeShortCircuit
    node
    R
    tspan_fault
end
NodeShortCircuit(;node,R,sc_timespan) = NodeShortCircuit(node,R,sc_timespan)

function (nsc::NodeShortCircuit)(powergrid)
    # Currently this assumes that the lines are PiModels...
    lines = copy(powergrid.lines)
    l_idx = findfirst(l -> (l.from == nsc.node || l.to == nsc.node) && l isa PiModelLine, lines)

    if l_idx == nothing
        @warn "Node needs to be connected to a PiModelLine to implement NodeShortCircuit"
        return nothing
    end

    if lines[l_idx].from == nsc.node
        lines[l_idx] = PiModelLine(;from=lines[l_idx].from, to=lines[l_idx].to, y = lines[l_idx].y, y_shunt_km = 1/nsc.R,  y_shunt_mk = lines[l_idx].y_shunt_mk)
    elseif lines[l_idx].to == nsc.node
        lines[l_idx] = PiModelLine(;from=lines[l_idx].from, to=lines[l_idx].to, y = lines[l_idx].y, y_shunt_km = lines[l_idx].y_shunt_km,  y_shunt_mk = 1/nsc.R)
    end
    PowerGrid(powergrid.nodes, lines)
end

"""
```Julia
simulate(nsc::NodeShortCircuit, powergrid, x1, timespan)
```
Simulates a [`NodeShortCircuit`](@ref)
"""
<<<<<<< HEAD

=======
>>>>>>> cleanup
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
        x3 = find_valid_initial_condition(nsc_powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(powergrid)
        integrator.u = x3
    end

    cb1 = DiscreteCallback(((u,t,integrator) -> t in nsc.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in nsc.tspan_fault[2]), regularState)
    sol = solve(problem, Rodas4(autodiff=false), callback = CallbackSet(cb1, cb2), tstops=[nsc.tspan_fault[1];nsc.tspan_fault[2]])
    return PowerGridSolution(sol, powergrid)
end
<<<<<<< HEAD
=======
export simulate
>>>>>>> cleanup
