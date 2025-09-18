"""
    line_between_slacks(l)

This function puts the line model between two slack nodes. The simulation consists of 3 phases:
 - 0-1s: Both slack nodes show the same voltage
 - 1-2s: Same magnitude, `dst` end +1 degree phase lead
 - 2-3s: Same magnitude, `dst` end -1 degree phase lag
 - 3-4s: `dst` end +10% magnitutde, same phase
 - 4-5s: `dst` end +10% magnitude, `dst` end +1 degree phase lead
 - 5-6s: `dst` end +10% magnitude, `dst` end -1 degree phase lag
"""
function line_between_slacks(edgef)
    src = compile_bus(SlackDifferential(name=:slack_src))
    dst = compile_bus(SlackDifferential(name=:slack_dst))
    g = path_graph(2)
    nw = Network(g, [src, dst], edgef)
    u0 = NWState(nw)
    any(isnan.(uflat(u0))) && error("Initial conditions contain NaNs")
    any(isnan.(pflat(u0))) && @warn "Parameters contain NaNs"

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

    plotspec = OrderedDict(
        "voltage angle" => OrderedDict(
            "angle at dst" => VIndex(2, :busbar₊u_arg)),
        "voltage magnitude" => OrderedDict(
            "magnitude at dst" => VIndex(2, :busbar₊u_mag)),
        "src end active power" => OrderedDict(
            "line injection toward src" => EIndex(1, :src₊P)),
        "dst end active power" => OrderedDict(
            "line injection toward dst" => EIndex(1, :dst₊P)),
        "src end reactive power" => OrderedDict(
            "line injection toward src" => EIndex(1, :src₊Q)),
        "dst end reactive power" => OrderedDict(
            "line injection toward dst" => EIndex(1, :dst₊Q))
    )

    return TrajectoriesOfInterest(sol, plotspec)
end

function bus_on_slack(busf; tmax=6, toilength=1000, argscale=1, magscale=1)
    slack = compile_bus(SlackDifferential(name=:slack_src))
    @named branch = PiLine(X=0.1)
    edgef = compile_line(MTKLine(branch))
    g = path_graph(2)

    nw = Network(g, [slack, busf], edgef)
    u0 = NWState(nw)
    if any(isnan.(uflat(u0)))
        # show(stderr, MIME"text/plain"(), u0)
        # println(stderr)
        error("Initial conditions contain NaNs")

    end
    if any(isnan.(pflat(u0)))
        # show(stderr, MIME"text/plain"(), u0.p)
        # println(stderr)
        @warn "Parameters contain NaNs"
        printstyled("Parameters contain NaNs:\n", color=:yellow)
        for (s, val) in zip(SII.parameter_symbols(u0), pflat(u0))
            if isnan(val)
                printstyled(" - $s = NaN\n", color=:yellow)
            end
        end
    end

    tstops = collect(range(0,tmax, length=7))[2:end-1]
    affect = function(integrator)
        u = NWState(integrator)
        if integrator.t == tstops[1]
           _set_voltage(u, 1, 1,  1/360*2π*argscale)
        elseif integrator.t == tstops[2]
           _set_voltage(u, 1, 1, -1/360*2π*argscale)
        elseif integrator.t == tstops[3]
           _set_voltage(u, 1, 1+0.1*magscale, 0)
        elseif integrator.t == tstops[4]
           _set_voltage(u, 1, 1+0.1*magscale,  1/360*2π*argscale)
        elseif integrator.t == tstops[5]
           _set_voltage(u, 1, 1+0.1*magscale, -1/360*2π*argscale)
        end
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback(tstops, affect)

    prob = ODEProblem(nw, uflat(u0), (0,tmax), pflat(u0); callback=cb)
    sol = solve(prob, Rodas5P())

    plotspec = OrderedDict(
        "active power" => OrderedDict("injection from bus" => VIndex(2, :busbar₊P)),
        "reactive power" => OrderedDict("injection from bus" => VIndex(2, :busbar₊Q)),
        "voltage angle" => OrderedDict(
            "angle at bus" => VIndex(2, :busbar₊u_arg),
            "angle at slack" => VIndex(1, :busbar₊u_arg)),
        "voltage magnitude" => OrderedDict(
            "magnitude at bus" => VIndex(2, :busbar₊u_mag),
            "magnitude at slack" => VIndex(1, :busbar₊u_mag)))
        # ,
        # "frequency" => OrderedDict("frequency at bus" => VIndex(2, :busbar₊ω)))

    return TrajectoriesOfInterest(sol, plotspec; length=toilength)
end

# using PowerDynamics.ModelChecks: _set_voltage
function _set_voltage(u::NWState, i, mag, arg)
    u.v[i,:busbar₊u_r] = mag * cos(arg)
    u.v[i,:busbar₊u_i] = mag * sin(arg)
end

function OpenIPSL_SMIB(_bus1)
    # copy constructor and set vidxs
    bus1 = VertexModel(_bus1, vidx=1, name=:GEN1)

    S_b = 100e6

    bus2 = let
        P_0=50e6
        Q_0=10e6
        v_0=0.9919935
        S_p_re = P_0
        S_p_im = Q_0
        @named constantLoad = PSSE_Load(; S_b, P_0, Q_0, v_0, S_p_re, S_p_im, characteristic=2)
        busmodel = MTKBus(constantLoad; name=:LOAD)
        bm = compile_bus(busmodel, pf=pfPQ(P=-P_0/S_b, Q=-Q_0/S_b), vidx=2)
    end

    bus3 = let
        slack = pfSlack(V=1; name=:GEN2)
        compile_bus(slack, vidx=3)
    end

    bus4 = let
        @named pwFault = ConstantYLoad(B=0, G=0)
        busmodel = MTKBus(pwFault; name=:FAULT)
        faultbus = compile_bus(busmodel, vidx=4)

        enable = ComponentAffect([], [:pwFault₊B, :pwFault₊G]) do u, p, ctx
            p[:pwFault₊B] = -1
            p[:pwFault₊G] = -1
        end
        disable = ComponentAffect([], [:pwFault₊B, :pwFault₊G]) do u, p, ctx
            p[:pwFault₊B] = 0
            p[:pwFault₊G] = 0
        end
        enable_cb = PresetTimeComponentCallback(2, enable)
        disable_cb = PresetTimeComponentCallback(2.15, disable)
        set_callback!(faultbus, (enable_cb, disable_cb))
        faultbus
    end

    bus5 = let
        busmodel = MTKBus(; name=:SHUNT)
        compile_bus(busmodel, vidx=5)
    end

    # line template
    pwLine = MTKLine(PiLine(; name=:PwLine))

    line = compile_line(pwLine; name=:pwLine,
        src=:GEN1, dst=:LOAD,
        PwLine₊X=0.2, PwLine₊R=0.001)
    line1 = compile_line(pwLine; name=:pwLine1,
        src=:LOAD, dst=:SHUNT,
        PwLine₊X=0.1, PwLine₊R=0.0005)
    line2 = compile_line(pwLine; name=:pwLine2,
        src=:SHUNT, dst=:GEN2,
        PwLine₊X=0.1, PwLine₊R=0.0005)
    line3 = compile_line(pwLine; name=:pwLine3,
        src=:LOAD, dst=:FAULT,
        PwLine₊X=0.1, PwLine₊R=0.0005)
    line4 = compile_line(pwLine; name=:pwLine4,
        src=:FAULT, dst=:GEN2,
        PwLine₊X=0.1, PwLine₊R=0.0005)

    buses = [bus1, bus2, bus3, bus4, bus5]
    lines = [line, line1, line2, line3, line4]
    nw = Network(buses, lines; warn_order=false)

    # pfnw = powerflow_model(nw)
    # pfs = solve_powerflow(pfnw)

    s0 = initialize_from_pf(nw, tol=1e-7, nwtol=1e-7)

    prob = ODEProblem(nw, uflat(s0), (0, 10), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())

    # plotspec = OrderedDict(
    #     "active power" => OrderedDict("injection from bus" => VIndex(1, :busbar₊P)),
    #     "reactive power" => OrderedDict("injection from bus" => VIndex(1, :busbar₊Q)),
    #     "voltage angle" => OrderedDict(
    #         "angle at LOAD" => VIndex(2, :busbar₊u_arg),
    #         "angle at GEN2" => VIndex(3, :busbar₊u_arg)),
    #     "voltage magnitude" => OrderedDict(
    #         "magnitude at LOAD" => VIndex(2, :busbar₊u_mag),
    #         "magnitude at GEN2" => VIndex(3, :busbar₊u_mag)))

    # prob = ODEProblem(nw, uflat(u0), (0, 10), pflat(u0))
    # sol = solve(prob, Rodas5P())

    # return TrajectoriesOfInterest(sol, plotspec)
end
