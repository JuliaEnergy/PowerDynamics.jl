# PSSE GENCLS Model - Port of OpenIPSL.Electrical.Machines.PSSE.GENCLS
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.
#
# Description: Classic generator model that can also represent an infinite bus

@mtkmodel PSSE_GENCLS begin
    @structural_parameters begin
        infinite_bus = false  # When true, H=0 behavior (infinite bus)
    end
    @components begin
        terminal = Terminal()
        # outputs
        δout = RealOutput() # rotor angle
        ωout = RealOutput() # rotor speed [pu]
        v_mag_out = RealOutput() # terminal voltage [pu]
        Pout = RealOutput() # active power [pu]
        Qout = RealOutput() # reactive power [pu]
    end
    @parameters begin
        # Machine parameters (using OpenIPSL naming)
        M_b, [description="Machine base power rating [MVA]"]
        H=0, [description="Inertia constant [s]"]
        D=0, [description="Damping coefficient [pu]"]
        R_a=0, [description="Armature resistance [pu]"]
        X_d=0.2, [description="d-axis transient reactance [pu]"]

        # System base
        S_b, [description="System power basis [MVA]"]
        V_b, [description="System voltage basis [kV]"]
        ω_b, [description="System base frequency [rad/s]"]

        # Initial conditions and setpoints
        P_0, [guess=0, description="Initial active power [MW]"]
        # Q_0, v_0 and angle_0 not actually used
        # Q_0, [guess=0, description="Initial reactive power [Mvar]"]
        # v_0, [guess=1, description="Initial voltage magnitude [pu]"]
        # angle_0, [guess=1, description="Initial voltage angle [rad]"]
        fn=50, [description="System frequency [Hz]"]
    end
    @variables begin
        # State variables
        δ(t), [guess=0, description="Rotor angle [rad]"]
        ω(t), [guess=0, description="Rotor speed deviation [pu]"]
        eq(t), [guess=1, description="Constant emf behind transient reactance [pu]"]

        # Algebraic variables
        vd(t), [description="d-axis voltage [pu]"]
        vq(t), [description="q-axis voltage [pu]"]
        id(t), [guess=0, description="d-axis current [pu]"]
        iq(t), [guess=0, description="q-axis current [pu]"]

        # OpenIPSL-style variables for equation compatibility
        pir(t), [description="Real part of terminal current (OpenIPSL convention) [pu]"]
        pii(t), [description="Imaginary part of terminal current (OpenIPSL convention) [pu]"]
        pvr(t), [description="Real part of terminal voltage (OpenIPSL convention) [pu]"]
        pvi(t), [description="Imaginary part of terminal voltage (OpenIPSL convention) [pu]"]

        # Observables
        V(t), [description="Bus voltage magnitude [pu]"]
        anglev(t), [description="Bus voltage angle [rad]"]
        P(t), [description="Active power [pu]"]
        Q(t), [description="Reactive power [pu]"]
    end
    begin
        # Base conversion factor
        CoB = M_b / S_b

        # Initial value calculations (from OpenIPSL)
        #=
        p0 = P_0 / M_b  # Initial active power (machine base)
        q0 = Q_0 / M_b  # Initial reactive power (machine base)
        vr0 = v_0 * cos(angle_0)
        vi0 = v_0 * sin(angle_0)
        ir0 = (p0 * vr0 + q0 * vi0) / (vr0^2 + vi0^2)
        ii0 = (p0 * vi0 - q0 * vr0) / (vr0^2 + vi0^2)
        delta0 = atan(vi0 + R_a*ii0 + X_d*ir0, vr0 + R_a*ir0 - X_d*ii0)
        vd0 = vr0 * cos(π/2 - delta0) - vi0 * sin(π/2 - delta0)
        vq0 = vr0 * sin(π/2 - delta0) + vi0 * cos(π/2 - delta0)
        id0 = ir0 * cos(π/2 - delta0) - ii0 * sin(π/2 - delta0)
        iq0 = ir0 * sin(π/2 - delta0) + ii0 * cos(π/2 - delta0)
        vf0 = vq0 + R_a*iq0 + X_d*id0
        =#
    end
    @equations begin
        # Conversion between PowerDynamics and OpenIPSL conventions
        # PowerDynamics: currents into terminal, OpenIPSL: currents out of terminal
        pir ~ -terminal.i_r
        pii ~ -terminal.i_i
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i

        # Swing equation - conditional behavior based on H
        # equivalent to abs(H) > C.eps in OpenIPSL
        # H=0 -> inifinite bus (no dynamics)
        Dt(δ) ~ ifelse(H>1e-10,
            ω * 2π * fn,
            0
        )
        Dt(ω) ~ ifelse(H>1e-10,
            (P_0/S_b - P - D*ω) / (2*H),
            0
        )

        # Classical model assumption: constant emf
        Dt(eq) ~ 0

        # Voltage equations (from OpenIPSL)
        vq ~ eq - R_a*iq - X_d*id  # q-axis voltage equation
        vd ~ X_d*iq - R_a*id       # d-axis voltage equation

        # Park's transformation (same as in OpenIPSL (and rest of PD for that matter))
        [pir, pii] .~ -CoB * [sin(δ)  cos(δ); -cos(δ)  sin(δ)] * [id, iq]
        [pvr, pvi] .~ [sin(δ)  cos(δ); -cos(δ)  sin(δ)] * [vd, vq]

        # Power injections (identical to OpenIPSL)
        -P ~ pvr * pir + pvi * pii
        -Q ~ pvi * pir - pvr * pii

        # Terminal voltage magnitude and angle (using OpenIPSL variables)
        V ~ sqrt(pvr^2 + pvi^2)
        anglev ~ atan(pvi, pvr)

        # Outputs
        δout.u ~ δ
        ωout.u ~ ω
        v_mag_out.u ~ V
        Pout.u ~ P
        Qout.u ~ Q
    end
end
