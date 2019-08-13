using Pkg
Pkg.activate(@__DIR__)
#Pkg.instantiate()

using Plots
using PowerDynamics # Currently requires master! use "]add PowerDynamics#master"
using NetworkDynamics
using DifferentialEquations

busses = [VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
    VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
    VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  PQAlgebraic(S=-1.0-im*0.5),
  PQAlgebraic(S=-1.0-im*0.5),
  PQAlgebraic(S=-1.0-im*0.5)]
  #RLC_Load(R,L,C),RLC_Load(R,L,C),RLC_Load(R,L,C)]


lines = [
  PiModelLine(;from=1, to=4, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=5, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=3, to=6, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=3, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=1, to=2, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.)]

pg = PowerGrid(busses, lines)

nsc = NodeShortCircuit(R = 0.1, node = 1, sc_timespan=(1.,5.))

#Find a stable operating point
ic_guess = ones(systemsize(pg)) + 0. * randn(systemsize(pg))
ic_guess[[3,6,9]] .= 0.
ic_guess = find_valid_ic(rhs(pg),ic_guess)
problem = ODEProblem(rhs(pg),ic_guess,(0.,200.))
sol = solve(problem, Rodas4(autodiff=false))
op_point = sol[end]

s1, s2, s3 = simulate(nsc, pg, op_point, (0.,30.));

plt=plot()
plot!(plt, s1, vars=[3,6,9])
plot!(plt, s2, vars=[3,6,9])
plot!(plt, s3, vars=[3,6,9])
plot!(plt, xlim=(0.,30.))
#
# pgs = PowerDynamics.PowerGridSolution(sol, pg)
# pgs(200., 1:3, :ω)

# Alternative Design:


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

export NodeShortCircuit2

nsc2 = NodeShortCircuit2(R = 0.1, node = 1, sc_timespan=(1.,5.))

sol = simulate(nsc2, pg, op_point, timespan=(0.,30.));

plot(sol)
