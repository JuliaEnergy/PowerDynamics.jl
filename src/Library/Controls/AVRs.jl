@mtkmodel AVRFixed begin
    @components begin
        vf = RealOutput()
    end
    @parameters begin
        vf_fixed, [guess=1, description="Fixed field voltage"]
    end
    @equations begin
        vf.u ~ vf_fixed
    end
end

@mtkmodel AVRTypeI begin
    @structural_parameters begin
       vref_input=false
       tmeas_lag=true
    end
    @components begin
        vh = RealInput() # generator terminal voltage
        if vref_input
            vref = RealInput() # reference voltage
        end
        vf = RealOutput(guess=1) # field voltage
    end
    @parameters begin
        Ka, [description="amplifier gain"]
        Ke, [description="field circuit integral deviation"]
        Kf, [description="stabilizer gain"]
        Ta, [description="amplifier time constant"]
        Tf, [description="stabilizer time constant"]
        Te, [description="field circuit time constant"]
        if tmeas_lag
            Tr, [description="measurement time constant"]
        end
        vr_min, [description="minimum regulator voltage"]
        vr_max, [description="maximum regulator voltage"]
        Ae, [description="1st ceiling coeff"]
        Be, [description="2st ceiling coeff"]
        if !vref_input
            vref, [guess=1, description="Terminal voltag reference [Machine PU]"]
        end
    end
    begin
        _vref = vref_input ? vref.u : vref
    end
    @variables begin
        vr1(t), [guess=1,description="intenral amplifier state TODO: implement AVR limiter"]
        vr2(t), [guess=0,description="internal feedback state"]
        vm(t), [guess=1, description="terminal voltage measurement (lagged)"]
        vfceil(t), [description="ceiled field voltage"]
    end
    @equations begin
        vfceil ~ Ae * exp(Be * abs(vf.u))
        Dt(vf.u) ~ -(vf.u  * (Ke+vfceil) - vr1)/Te
        Dt(vr1) ~ (Ka*(_vref - vm - vr2 - Kf/Tf*vf.u) - vr1)/Ta
        Dt(vr2) ~ -(Kf/Tf*vf.u + vr2)/Tf
        if tmeas_lag
            Dt(vm) ~ (vh.u - vm)/Tr
        else
            vm ~ vh.u
        end
    end
end

function solve_ceilf(pair1, pair2)
    p = [pair1.first, pair1.second, pair2.first, pair2.second]
    f = (du, u, p) -> begin
        A,B = u
        v1, S1, v2, S2 = p
        du[1] = S1 - A*exp(B*v1)
        du[2] = S2 - A*exp(B*v2)
    end
    prob = NonlinearProblem(f, [0.01, 1], p)
    sol = solve(prob)
    if !SciMLBase.successful_retcode(sol.retcode)
        error("Did not finde solution for Ae and Be: retcode $(sol.retcode)")
    end
    (; Ae=sol[1], Be=sol[2])
end
