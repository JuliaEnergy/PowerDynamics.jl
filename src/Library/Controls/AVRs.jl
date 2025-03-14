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
       ceiling_function=:exponential
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
        E1, [description="1st ceiling voltage"]
        E2, [description="2nd ceiling voltage"]
        Se1, [description="1st ceiling saturation"]
        Se2, [description="2nd ceiling saturation"]
        if ceiling_function == :exponential
            Ae=_solve_Ae(E1=>Se1, E2=>Se2), [description="1st ceiling coeff"]
            Be=_solve_Be(E1=>Se1, E2=>Se2), [description="1st ceiling coeff"]
        end
        if !vref_input
            vref, [guess=1, description="Terminal voltag reference [Machine PU]"]
        end
    end
    begin
        _vref = vref_input ? vref.u : vref
        if ceiling_function ∉ (:exponential, :quadratic)
            error("Unknown ceiling function: $ceiling_function")
        end
    end
    @variables begin
        # we add an explicit state for vfout to set bounds
        vfout(t), [guess=1, bounds=(0,Inf), description="field voltage output"]
        vr(t), [guess=0, description="regulator voltage"]
        vm(t), [guess=1, description="terminal voltage measurement (lagged)"]
        vfceil(t), [description="ceiled field voltage"]
        amp_in(t), [description="amplifier input"]
        v_fb(t), [guess=0, description="feedback voltage"]
    end
    @equations begin
        # implementation after block diagram in milano
        if ceiling_function == :exponential
            vfceil ~ vfout * Ae * exp(Be * abs(vfout))
        elseif ceiling_function == :quadratic
            vfceil ~ quadratic_ceiling(abs(vfout), E1, E2, Se1, Se2)
        end
        if tmeas_lag
            Tr * Dt(vm) ~ vh.u - vm
        else
            vm ~ vh.u
        end

        Tf*Dt(v_fb) ~ Kf*Dt(vfout) - v_fb
        amp_in ~ Ka*(_vref - vm - v_fb)

        Ta*Dt(vr) ~ ifelse(
            ((vr > vr_max) & (amp_in > vr)) | ((vr < vr_min) & (amp_in < vr)),
            0,
            amp_in - vr)

        Te*Dt(vfout) ~ vr - vfceil - Ke*vfout

        # output
        vf.u ~ vfout
    end
end

function solve_ceilf(pair1, pair2; u0=[0.01, 1])
    p = [pair1.first, pair1.second, pair2.first, pair2.second]
    f = (du, u, p) -> begin
        A,B = u
        v1, S1, v2, S2 = p
        du[1] = S1 - A*exp(B*v1)
        du[2] = S2 - A*exp(B*v2)
    end
    prob = NonlinearProblem(f, u0, p)
    sol = solve(prob; verbose=false)
    if !SciMLBase.successful_retcode(sol.retcode)
        error("Did not finde solution for Ae and Be: retcode $(sol.retcode)")
    end
    (; Ae=sol[1], Be=sol[2])
end
_solve_Ae(pair1, pair2) = solve_ceilf(pair1, pair2)[1]
_solve_Be(pair1, pair2) = solve_ceilf(pair1, pair2)[2]


function quadratic_ceiling(x, E1, E2, Se1, Se2)
    # sq = sqrt(Se1/Se2)
    # Asq = (E1 - E2 * sq) / (1 - sq)
    # Bsq = Se2 /(E2 - Asq)^2

    # below is the definition from RMSPowerSims
    sq = sqrt((E1 * Se1) / (E2 * Se2))
    Asq = (E1 - E2 * sq) / (1 - sq)
    Bsq = (E2 * Se2) / ((E2 - Asq)^2)

    ifelse(x > Asq, Bsq * (x - Asq)^2, 0.0)
end


@mtkmodel AVRTypeIS begin
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
        A, [description="saturation constant Ag"]
        B, [description="saturation constant Bg"]
        if !vref_input
            vref, [guess=1, description="Terminal voltag reference [Machine PU]"]
        end
    end
    begin
        _vref = vref_input ? vref.u : vref
    end
    @variables begin
        # we add an explicit state for vfout to set bounds
        vfout(t), [guess=1, bounds=(0,Inf), description="field voltage output"]
        vr(t), [guess=0, description="regulator voltage"]
        vm(t), [guess=1, description="terminal voltage measurement (lagged)"]
        vfceil(t), [description="ceiled field voltage"]
        amp_in(t), [description="amplifier input"]
        v_fb(t), [guess=0, description="feedback voltage"]
    end
    @equations begin
        # implementation after block diagram in milano / PowerFactory
        vfceil ~ A * exp(B * vfout) * vfout
        if tmeas_lag
            Tr * Dt(vm) ~ vh.u - vm
        else
            vm ~ vh.u
        end

        Tf*Dt(v_fb) ~ Kf*Dt(vfout) - v_fb
        amp_in ~ Ka*(_vref - vm - v_fb)

        Ta*Dt(vr) ~ ifelse(
            ((vr > vr_max) & (amp_in > vr)) | ((vr < vr_min) & (amp_in < vr)),
            0,
            amp_in - vr)

        Te*Dt(vfout) ~ vr - vfceil - Ke*vf.u

        # output
        vf.u ~ vfout
    end
end