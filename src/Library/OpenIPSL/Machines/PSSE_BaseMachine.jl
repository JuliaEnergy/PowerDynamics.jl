# PSSE Basemachine Model - Port of OpenIPSL.Electrical.Machines.PSSE.baseMachine
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.
#

@mtkmodel PSSE_BaseMachine begin
    @components begin
        terminal = Terminal()

        # Output interfaces
        SPEED_out = RealOutput()  # Machine speed deviation from nominal [pu]
        ETERM_out = RealOutput()  # Machine terminal voltage [pu]
        PELEC_out = RealOutput()  # Machine electrical power
        ISORCE_out = RealOutput() # Machine source current [pu]
        ANGLE_out = RealOutput()  # Machine relative rotor angle
        XADIFD_out = RealOutput() # Machine field current [pu]
        # initial state outputs not necessary!
        # PMECH0_out = RealOutput() # Initial value of machine electrical power
        # EFD0_out = RealOutput()   # Initial generator main field voltage [pu]
    end

    @parameters begin
        # Machine parameters (free parameters - need guess values)
        M_b, [description="Machine base power [MVA]"]
        Tpd0, [description="d-axis transient open-circuit time constant [s]"]
        Tppd0, [description="d-axis sub-transient open-circuit time constant [s]"]
        Tppq0, [description="q-axis sub-transient open-circuit time constant [s]"]
        H, [description="Inertia constant [s]"]
        D, [description="Speed damping [pu]"]
        Xd, [description="d-axis reactance [pu]"]
        Xq, [description="q-axis reactance [pu]"]
        Xpd, [description="d-axis transient reactance [pu]"]
        Xppd, [description="d-axis sub-transient reactance [pu]"]
        Xppq, [description="q-axis sub-transient reactance [pu]"]
        Xl, [description="leakage reactance [pu]"]
        S10, [description="Saturation factor at 1.0 pu [pu]"]
        S12, [description="Saturation factor at 1.2 pu [pu]"]

        # System parameters (free parameters from pfComponent)
        S_b, [description="System power basis [MVA]"]
        fn=50, [description="System frequency [Hz]"]

        # Fixed parameters (with default values)
        R_a=0, [description="Armature resistance [pu]"]
        # w0=0, [description="Initial speed deviation from nominal [pu]"]
    end

    @variables begin
        # State variables
        w(t), [guess=0, description="Machine speed deviation [pu]"]
        delta(t), [guess=0, description="Rotor angle [rad]"]

        # Algebraic variables
        Vt(t), [guess=1.0, description="Bus voltage magnitude [pu]"]
        anglev(t), [guess=0, description="Bus voltage angle [rad]"]
        I(t), [guess=1.0, description="Terminal current magnitude [pu]"]
        anglei(t), [guess=0, description="Terminal current angle [rad]"]
        P(t), [guess=0.4, description="Active power (system base) [pu]"]
        Q(t), [guess=0.05, description="Reactive power (system base) [pu]"]
        Te(t), [guess=0.4, description="Electrical torque [pu]"]
        id(t), [guess=0, description="d-axis armature current [pu]"]
        iq(t), [guess=1.0, description="q-axis armature current [pu]"]
        ud(t), [guess=0, description="d-axis terminal voltage [pu]"]
        uq(t), [guess=1.0, description="q-axis terminal voltage [pu]"]

        # OpenIPSL-style variables for equation compatibility
        pir(t), [guess=0, description="Real part of terminal current (OpenIPSL convention) [pu]"]
        pii(t), [guess=0, description="Imaginary part of terminal current (OpenIPSL convention) [pu]"]
        pvr(t), [guess=1, description="Real part of terminal voltage (OpenIPSL convention) [pu]"]
        pvi(t), [guess=0, description="Imaginary part of terminal voltage (OpenIPSL convention) [pu]"]

        # Input/parameter variables
        pmech(t), [description="mechanical power [pu]"]
        efd(t), [description="field voltage [pu]"]
    end

    begin
        # Derived parameters (used in equations)
        CoB = M_b/S_b # Base conversion factor
    end

    @equations begin
        # Conversion between PowerDynamics and OpenIPSL conventions
        pir ~ -terminal.i_r
        pii ~ -terminal.i_i
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i

        # Interfacing outputs with internal variables (from OpenIPSL line 109-113)
        ANGLE_out.u ~ delta
        SPEED_out.u ~ w
        ETERM_out.u ~ Vt
        PELEC_out.u ~ P/CoB

        # Current and voltage transformations (from OpenIPSL line 114-115)
        [pir, pii] .~ -CoB*[sin(delta)  cos(delta); -cos(delta)  sin(delta)] * [id, iq]
        [pvr, pvi] .~ [sin(delta)  cos(delta); -cos(delta)  sin(delta)] * [ud, uq]

        # Power calculations (from OpenIPSL line 116-117)
        -P ~ pvr*pir + pvi*pii
        -Q ~ pvi*pir - pvr*pii

        # Terminal voltage and current magnitudes/angles (from OpenIPSL line 118-121)
        Vt ~ sqrt(pvr^2 + pvi^2)
        anglev ~ atan(pvi, pvr)
        I ~ sqrt(pii^2 + pir^2)
        anglei ~ atan(pii, pir)

        # Swing equations (from OpenIPSL line 122-123)
        Dt(w) ~ ((pmech - D*w)/(w + 1) - Te)/(2*H)
        Dt(delta) ~ 2*π*fn*w
    end
end
