"""
    AVRFixed

Trivial AVR that holds the field voltage at a fixed parameter value.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
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

"""
    AVRTypeI

IEEE Type I excitation system with amplifier, stabilizer feedback, field circuit, and ceiling function.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel AVRTypeI begin
    @structural_parameters begin
       vref_input=false
       tmeas_lag=true
       anti_windup = true
       ceiling_function=:exponential
    end
    @components begin
        v_mag = RealInput() # generator terminal voltage
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
        vr1(t), [guess=0, description="regulator voltage before Limiter"]
        v_fb(t), [guess=0, description="feedback voltage"]
    end
    @equations begin
        # implementation after block diagram in milano
        if ceiling_function == :exponential
            vfceil ~ vfout * EXP_SE(abs(vfout), Se1, Se2, E1, E2)
        elseif ceiling_function == :quadratic
            vfceil ~ abs(vfout) * QUAD_SE(abs(vfout), Se1, Se2, E1, E2)
        end
        if tmeas_lag
            Tr * Dt(vm) ~ v_mag.u - vm
        else
            vm ~ v_mag.u
        end

        Tf*Dt(v_fb) ~ Kf*Dt(vfout) - v_fb

        if anti_windup
            amp_in ~ Ka*(_vref - vm - v_fb)
            Ta*Dt(vr) ~ ifelse(
                ((vr > vr_max) & (amp_in > vr)) | ((vr < vr_min) & (amp_in < vr)),
                0,
                amp_in - vr)
        else
            Ta*Dt(vr1) ~ Ka*(_vref - vm - v_fb) - vr1
            vr ~ max(vr_min, min(vr_max, vr1))
        end

        Te*Dt(vfout) ~ vr - vfceil - Ke*vfout

        # output
        vf.u ~ vfout
    end
end


#=
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
        vr1(t), [guess=0, description="regulator voltage before Limiter"]
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
        #amp_in ~ Ka*(_vref - vm - v_fb)

        #Ta*Dt(vr) ~ ifelse(
           # ((vr > vr_max) & (amp_in > vr)) | ((vr < vr_min) & (amp_in < vr)),
           # 0,
            #amp_in - vr)
        Ta*Dt(vr1) ~ Ka*(_vref - vm - v_fb) - vr1
        vr ~ max(vr_min, min(vr_max, vr1))

        Te*Dt(vfout) ~ vr - vfceil - Ke*vf.u

        # output
        vf.u ~ vfout
    end
end
=#
