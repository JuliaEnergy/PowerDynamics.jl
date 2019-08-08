export NodeShortCircuit

using NetworkDynamics
using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

@Base.kwdef struct NodeShortCircuit
    R
    node
    sc_timespan
end

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


function simulate(nsc::NodeShortCircuit, powergrid, x1, timespan)
    sc_timespan = nsc.sc_timespan
    @assert timespan[1] < sc_timespan[1]
    @assert timespan[2] > sc_timespan[2]
    nsc_powergrid = nsc(powergrid)

    # Integrate to fault
    prob1 = ODEProblem(rhs(powergrid), x1, (timespan[1], sc_timespan[1]))
    sol1 = solve(prob1, Rodas4(autodiff=false))

    # Integrate the fault state
    x2 = find_valid_ic(rhs(nsc_powergrid), sol1[end]) # Jump the state to be valid for the new system.
    prob2 = ODEProblem(rhs(nsc_powergrid), x2, sc_timespan)
    sol2 = solve(prob2, Rodas4(autodiff=false))

    # Integrate after fault
    x3 = find_valid_ic(rhs(powergrid), sol2[end]) # Jump the state to be valid for the new system.
    prob3 = ODEProblem(rhs(powergrid), x3, (sc_timespan[2], timespan[2]))
    sol3 = solve(prob3, Rodas4(autodiff=false))

    sol1, sol2, sol3
end
