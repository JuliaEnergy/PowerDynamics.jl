function OpenIPSL_SMIB(_bus1)
    # copy constructor and set vidxs
    bus1 = VertexModel(_bus1, vidx=1, name=:GEN1)

    S_b = 100e6

    bus2 = let
        # pf results for setting up pf model
        P_0=50e6
        Q_0=10e6
        v_0=0.9919935

        @named constantLoad = PSSE_Load(; v_0, S_b, characteristic=2)
        busmodel = MTKBus(constantLoad; name=:LOAD)
        bm = compile_bus(busmodel, pf=pfPQ(P=-P_0/S_b, Q=-Q_0/S_b), vidx=2)
    end

    bus3 = let
        # OpenIPSL infinite bus parameters from SMIB base class
        H = 0                    # H=0 makes it behave like infinite bus
        M_b = 100e6
        X_d = 0.2               # Internal impedance
        D = 0
        V_b = 400e3
        ω_b = 2π*50

        # pf results, just used for pf modek
        # P_0 = 10.017110e6       # From OpenIPSL SMIB.mo
        # Q_0 = 8.006544e6        # From OpenIPSL SMIB.mo
        v_0 = 1.0
        angle_0 = 0.0           # From OpenIPSL SMIB.mo

        @named gencls_inf = PSSE_GENCLS(; S_b, V_b, ω_b, H, M_b, X_d, D)
        busmodel = MTKBus(gencls_inf; name=:GEN2)
        compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0), vidx=3)
    end

    bus4 = let
        @named pwFault = ConstantYLoad(B=0, G=0, allow_zero_conductance=true)
        busmodel = MTKBus(pwFault; name=:FAULT)
        faultbus = compile_bus(busmodel, vidx=4)

        enable = ComponentAffect([], [:pwFault₊B, :pwFault₊G]) do u, p, ctx
            p[:pwFault₊B] = 1e10
            p[:pwFault₊G] = 1e10
            # p[:pwFault₊B] = 10
            # p[:pwFault₊G] = 10
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

    s0 = initialize_from_pf(nw, subverbose=true)

    prob = ODEProblem(nw, uflat(s0), (0, 10), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end

function ref_rms_error(sol, csv, idx, col)
    t = csv[:, "time"]
    _ref = csv[:, col]

    # we find alle t in ti, where the jump occurs
    # at these points, we ignore the error because its hard to define with left/right limit
    Δtmean = mean(diff(t))
    _ti_jump = findall(Δ -> Δ < Δtmean/10000, diff(t))
    ti_jump = unique!(sort!(vcat(_ti_jump, _ti_jump .+ 1)))

    deleteat!(t, ti_jump)
    deleteat!(_ref, ti_jump)

    _sim = sol(t, idxs=idx).u

    norm(_ref .- _sim) / sqrt(length(_ref))
end
