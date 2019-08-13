export NodeShortCircuit2

mutable struct pg_switch
    original_pg::PowerGrid
    fault_pg::PowerGrid
    fault::Bool
end

function pg_switch(pg::PowerGrid, fault)
    pg_switch(pg, fault(pg), false)
end

function (pgs::pg_switch)(u, t, integrator)
    println("time $t") # To make sure the call back is really called
    if ! pgs.fault
        integrator.f = rhs(pgs.fault_pg)
        integrator.u = find_valid_ic(integrator.f, integrator.u)
        # without the next lines I get wrong results, with them I get an
        # infinite loop as the callback is triggered over and over.
        # reinit!(integrator, t0=t)
        # step!(integrator)
        pgs.fault = true
    else
        integrator.f = rhs(pgs.original_pg)
        integrator.u = find_valid_ic(integrator.f, integrator.u)
        # without the next lines I get wrong results, with them I get an
        # infinite loop as the callback is triggered over and over.
        # reinit!(integrator, t0=t)
        # step!(integrator)
        pgs.fault = false
    end
end


@Base.kwdef struct NodeShortCircuit2
    R
    node
    sc_timespan
end

function (nsc::NodeShortCircuit2)(powergrid)
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

import PowerDynamics.simulate

function simulate(nsc::NodeShortCircuit2, powergrid, x0; timespan)
    @assert timespan[1] < nsc.sc_timespan[1]
    @assert timespan[2] > nsc.sc_timespan[2]
    pgs = pg_switch(powergrid, nsc)

    cb = FunctionCallingCallback(pgs; funcat=nsc.sc_timespan, func_start = false)

    prob = ODEProblem(rhs(powergrid), x0, timespan, callback=cb)

    solve(prob, Rodas4(autodiff=false))
end
