module ModelChecks

using Graphs: SimpleGraph, add_edge!, path_graph
using NetworkDynamics: Network, NWState, uflat, pflat, vidxs, eidxs
using OrdinaryDiffEq: ODEProblem, solve, Rodas5P, auto_dt_reset!
using DiffEqCallbacks: PresetTimeCallback
using ModelingToolkit: @named
using Makie: lines, Figure, Axis, axislegend, lines!

using ..OpPoDyn
using ..OpPoDyn.Library

export line_between_slacks, bus_on_slack

"""
    line_between_slacks(l::Line)

This function puts the line model between two slack nodes. The simulation consists of 3 phases:
 - 0-1s: Both slack nodes show the same voltage
 - 1-2s: Same magnitude, `dst` end +1 degree phase lead
 - 2-3s: Same magnitude, `dst` end -1 degree phase lag
 - 3-4s: `dst` end +10% magnitutde, same phase
 - 4-5s: `dst` end +10% magnitude, `dst` end +1 degree phase lead
 - 5-6s: `dst` end +10% magnitude, `dst` end -1 degree phase lag
"""
function line_between_slacks(line::Line)
    src = Bus(SlackDifferential(name=:slack_src)).compf
    dst = Bus(SlackDifferential(name=:slack_dst)).compf
    edgef = line.compf
    g = path_graph(2)
    nw = Network(g, [src, dst], edgef)
    u0 = NWState(nw)
    any(isnan.(uflat(u0))) && error("Initial conditions contain NaNs")
    any(isnan.(pflat(u0))) && error("Parameters contain NaNs")

    affect = function(integrator)
        u = NWState(integrator)
        if integrator.t == 1
           _set_voltage(u, 2, 1,  1/360*2π)
        elseif integrator.t == 2
           _set_voltage(u, 2, 1, -1/360*2π)
        elseif integrator.t == 3
           _set_voltage(u, 2, 1.1, 0)
        elseif integrator.t == 4
           _set_voltage(u, 2, 1.1,  1/360*2π)
        elseif integrator.t == 5
           _set_voltage(u, 2, 1.1, -1/360*2π)
        end
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback([1,2,3,4,5], affect)

    prob = ODEProblem(nw, uflat(u0), (0,6), pflat(u0); callback=cb)
    sol = solve(prob, Rodas5P())

    ts = range(0, 6, length=1000)
    P1 = sol(ts, idxs=vidxs(nw, 1, :busbar₊P))
    Q1 = sol(ts, idxs=vidxs(nw, 1, :busbar₊Q))
    P2 = sol(ts, idxs=vidxs(nw, 2, :busbar₊P))
    Q2 = sol(ts, idxs=vidxs(nw, 2, :busbar₊Q))
    srcP = sol(ts, idxs=eidxs(nw, 1, :src₊P))
    dstP = sol(ts, idxs=eidxs(nw, 1, :dst₊P))
    srcQ = sol(ts, idxs=eidxs(nw, 1, :src₊Q))
    dstQ = sol(ts, idxs=eidxs(nw, 1, :dst₊Q))

    fig = Figure(resolution=(2000,1000))
    ax1 = Axis(fig[1,1], title="src end active power")
    lines!(ax1, P1.t, P1.u; label="injection at bus 1")
    lines!(ax1, srcP.t, srcP.u; label="line injection toward bus")
    axislegend(ax1)

    ax2 = Axis(fig[1,2], title="dst end active power")
    lines!(ax2, P2.t, P2.u; label="injection at bus 2")
    lines!(ax2, dstP.t, dstP.u; label="line injection toward bus")
    axislegend(ax2)

    ax3 = Axis(fig[2,1], title="src end reactive power")
    lines!(ax3, Q1.t, Q1.u; label="injection at bus 1")
    lines!(ax3, srcQ.t, srcQ.u; label="line injection toward bus")
    axislegend(ax3)

    ax4 = Axis(fig[2,2], title="dst end reactive power")
    lines!(ax4, Q2.t, Q2.u; label="injection at bus 2")
    lines!(ax4, dstQ.t, dstQ.u; label="line injection toward bus")
    axislegend(ax4)
    fig
end

function bus_on_slack(bus::Bus)
    slack = Bus(SlackDifferential(name=:slack_src)).compf
    @named branch = DynawoPiLine(XPu=0.04189)
    edgef = Line(LineModel(branch)).compf
    busf = bus.compf
    g = path_graph(2)

    nw = Network(g, [slack, busf], edgef)
    u0 = NWState(nw)
    any(isnan.(uflat(u0))) && error("Initial conditions contain NaNs")
    any(isnan.(pflat(u0))) && error("Parameters contain NaNs")

    affect = function(integrator)
        u = NWState(integrator)
        if integrator.t == 1
           _set_voltage(u, 1, 1,  1/360*2π)
        elseif integrator.t == 2
           _set_voltage(u, 1, 1, -1/360*2π)
        elseif integrator.t == 3
           _set_voltage(u, 1, 1.1, 0)
        elseif integrator.t == 4
           _set_voltage(u, 1, 1.1,  1/360*2π)
        elseif integrator.t == 5
           _set_voltage(u, 1, 1.1, -1/360*2π)
        end
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback([1,2,3,4,5], affect)

    prob = ODEProblem(nw, uflat(u0), (0,6), pflat(u0); callback=cb)
    sol = solve(prob, Rodas5P())

    ts = range(0, 6, length=1000)
    P = sol(ts, idxs=vidxs(nw, 2, :busbar₊P))
    Q = sol(ts, idxs=vidxs(nw, 2, :busbar₊Q))
    δ = sol(ts, idxs=vidxs(nw, 2, :busbar₊u_arg))
    δslack = sol(ts, idxs=vidxs(nw, 1, :busbar₊u_arg))

    fig = Figure(size=(2000,1000))
    ax1 = Axis(fig[1,1], title="active power")
    lines!(ax1, P.t, P.u; label="injection at bus")
    axislegend(ax1)

    ax2 = Axis(fig[1,2], title="voltage angle")
    lines!(ax2, δ.t, δ.u; label="angle at bus")
    # lines!(ax2, δslack.t, δslack.u; label="angle at slack")
    axislegend(ax2)
    fig
end


# using OpPoDyn.ModelChecks: _set_voltage
function _set_voltage(u::NWState, i, mag, arg)
    u.v[i,:busbar₊u_r] = mag * cos(arg)
    u.v[i,:busbar₊u_i] = mag * sin(arg)
end

end # module
