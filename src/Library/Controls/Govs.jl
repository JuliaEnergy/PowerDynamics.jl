@mtkmodel GovFixed begin
    @components begin
        τ_m = RealOutput()
    end
    @parameters begin
        τ_m_fixed, [guess=0, description="Fixed mechanical torque"]
    end
    @equations begin
        τ_m.u ~ τ_m_fixed
    end
end

_clamp(u, u_min, u_max) = max(min(u, u_max), u_min)

# from Milano P. 359
@mtkmodel TurbineGovTypeI begin
    @structural_parameters begin
        ω_ref_input=false
        p_ref_input=false
    end
    @components begin
        ω_meas = RealInput()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if p_ref_input
            p_ref = RealInput()
        end
        τ_m = RealOutput()
    end
    @parameters begin
       if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
       end
       if !p_ref_input
            p_ref, [guess=1, description="Reference power [Machine PU]"]
       end
       p_min, [description="Minimum turbine output [Machine PU]"]
       p_max, [description="Maximum turbine output [Machine PU]"]
       R, [description="Govenor droop [Machine PU]"]
       # TODO: check Tc servo Ts governor
       Tc, [description="Servo time constant [s]"]
       Ts, [description="Govenor time constant [s]"]
       T3, [description="Transient time constant 3 [s]"]
       T4, [description="Transient time constant 4 [s]"]
       T5, [description="Transient time constant 5 [s]"]
    end
    @variables begin
        p_droop(t), [description="P after droop (not limited)"]
        p_lim(t), [description="limited p"]
        xg1(t)
        xg2(t)
        xg3(t)
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _p_ref = p_ref_input ? p_ref.u : p_ref
    end
    @equations begin
        p_droop ~ p_ref + 1/R * (_ω_ref - ω_meas)
        p_lim ~ _clamp(p_droop, p_min, p_max)
        Ts * Dt(xg1) ~ p_lim - xg1
        Tc * Dt(xg2) ~ (1-T3/Tc)*xg1 - xg2
        T5 * Dt(xg3) ~ (1-T4/T5)*(xg2 + T3/Tc*xg1) - xg3
        Dt(τ_m.u) ~ xg3 + T4/T5*(xg2 + T3/Tc*xg1)
    end
end

# from PSD.jl https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/turbine_gov/#TGOV1-[SteamTurbineGov1]
@mtkmodel TGOV1 begin
    @structural_parameters begin
        ω_ref_input=false
        p_ref_input=false
    end
    @components begin
        ω_meas = RealInput()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if p_ref_input
            p_ref = RealInput()
        end
        τ_m = RealOutput()
    end
    @parameters begin
       if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
       end
       if !p_ref_input
            p_ref, [guess=1, description="Reference power [Machine PU]"]
       end
       V_min, [description="Valve min position"]
       V_max, [description="Valve max position"]
       R, [description="Govenor droop [Machine PU]"]
       T1, [description="Transient time constant 1 [s]"]
       T2, [description="Transient time constant 2 [s]"]
       T3, [description="Transient time constant 3 [s]"]
       DT, [description="Turbine Damping"]
    end
    @variables begin
        ref_sig(t), [description="Internal reference signal"]
        xg1(t), [guess=0]
        xg1_sat(t)
        xg2(t), [guess=0]
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _p_ref = p_ref_input ? p_ref.u : p_ref
    end
    @equations begin
        ref_sig ~ 1/R*(_p_ref  - (ω_meas.u - _ω_ref))
        T1 * Dt(xg1) ~ (ref_sig - xg1)
        xg1_sat ~ _clamp(xg1, V_min, V_max)
        T3 * Dt(xg2) ~ xg1_sat*(1-T2/T3) - xg2
        # TODO: check units, might need multiplication by omega ref to get from p to tau
        τ_m.u ~ xg2 + T2/T3*xg1_sat - DT*(ω_meas.u - _ω_ref)
    end
end
