# PSSE Load Model - Port of OpenIPSL.Electrical.Loads.PSSE.Load
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.
#
# Description: PSSE Load model with ZIP characteristics, voltage dependency, and load variation capability
# Reference: PSS/E Manual

@mtkmodel PSSE_Load begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        # Power flow initial conditions (from pfComponent)
        P_0, [description="Initial active power [MW]"]
        Q_0, [description="Initial reactive power [Mvar]"]
        v_0, [description="Initial voltage magnitude [pu]"]
        # angle_0, [description="Initial voltage angle [rad]"]
        S_b, [description="System base power [MVA]"]

        # ZIP Load Components (from OpenIPSL Load model)
        S_p_re, [guess=0, description="Real part of original constant power load [MW]"]
        S_p_im, [guess=0, description="Imaginary part of original constant power load [Mvar]"]
        S_i_re=0.0, [description="Real part of original constant current load [MW]"]
        S_i_im=0.0, [description="Imaginary part of original constant current load [Mvar]"]
        S_y_re=0.0, [description="Real part of original constant admittance load [MW]"]
        S_y_im=0.0, [description="Imaginary part of original constant admittance load [Mvar]"]

        # Load Transfer Fractions
        a_re=1.0, [description="Real part of load transfer fraction for constant current load"]
        a_im=0.0, [description="Imaginary part of load transfer fraction for constant current load"]
        b_re=0.0, [description="Real part of load transfer fraction for constant admittance load"]
        b_im=1.0, [description="Imaginary part of load transfer fraction for constant admittance load"]

        # Voltage Dependency Characteristics
        PQBRAK=0.7, [description="Constant power characteristic threshold [pu]"]
        characteristic=1, [description="Voltage dependency characteristic type (1 or 2)"]

        # Load Variation
        d_P=0.0, [description="Active load variation [pu] (default 0 for no variation). Reactive load variation d_Q is automatically calculated to maintain constant power factor."]
    end
    @variables begin
        # OpenIPSL-style variables for equation compatibility
        pir(t), [guess=0, description="Real part of terminal current (OpenIPSL convention) [pu]"]
        pii(t), [guess=0, description="Imaginary part of terminal current (OpenIPSL convention) [pu]"]
        pvr(t), [guess=1, description="Real part of terminal voltage (OpenIPSL convention) [pu]"]
        pvi(t), [guess=0, description="Imaginary part of terminal voltage (OpenIPSL convention) [pu]"]

        # State variables
        v(t), [guess=1, description="Bus voltage magnitude [pu]"]
        angle(t), [guess=0, description="Bus voltage angle [rad]"]
        P(t), [guess=0, description="Active power consumption [pu]"]
        Q(t), [guess=0, description="Reactive power consumption [pu]"]

        # Voltage dependency factors
        kP(t), [guess=1, description="Constant power load factor"]
        kI(t), [guess=1, description="Constant current load factor"]
    end
    begin
        # Complex ZIP load components (from baseLoad protected parameters)
        # Convert MW/Mvar to pu and apply load transfer fractions
        S_P_re = ((1 - a_re - b_re) * S_p_re) / S_b
        S_P_im = ((1 - a_im - b_im) * S_p_im) / S_b
        S_I_re = (S_i_re + (a_re * S_p_re / v_0)) / S_b
        S_I_im = (S_i_im + (a_im * S_p_im / v_0)) / S_b
        S_Y_re = (S_y_re + (b_re * S_p_re / v_0^2)) / S_b
        S_Y_im = (S_y_im + (b_im * S_p_im / v_0^2)) / S_b

        # Voltage dependency coefficients for characteristic 2
        a2 = 1.502  # Constant current load coefficient
        b2 = 1.769  # Constant current load exponent
        a0 = 0.4881 # Constant power load base coefficient
        a1 = -0.4999 # Constant power load cosine coefficient
        b1 = 0.1389  # Constant power load sine coefficient
        wp = 3.964   # Constant power load frequency

        # Initial conditions (from baseLoad)
        p0 = (S_I_re * v_0 + S_Y_re * v_0^2 + S_P_re)
        q0 = (S_I_im * v_0 + S_Y_im * v_0^2 + S_P_im)
        # vr0 = v_0 * cos(angle_0)
        # vi0 = v_0 * sin(angle_0)
        # ir0 = (p0 * vr0 + q0 * vi0) / (vr0^2 + vi0^2)
        # ii0 = (p0 * vi0 - q0 * vr0) / (vr0^2 + vi0^2)

        # Load variation (maintains constant power factor like OpenIPSL)
        PF = ifelse(abs(q0) <= 1e-10, 1.0, p0/q0)  # Power factor ratio
        d_Q = (p0 + d_P)/PF - q0  # Derived reactive load variation
    end
    @equations begin
        # Conversion between PowerDynamics and OpenIPSL conventions
        pir ~ -terminal.i_r
        pii ~ -terminal.i_i
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i

        # Voltage magnitude and angle calculations
        v ~ sqrt(pvr^2 + pvi^2)
        angle ~ atan(pvi, pvr)

        # Power calculations
        P ~ pvr * pir + pvi * pii
        Q ~ (-pvr * pii) + pvi * pir

        # Voltage dependency characteristics (from baseLoad)
        # Characteristic 1: Parabolic behavior below PQBRAK
        kP ~ ifelse(characteristic == 1,
            ifelse(v < PQBRAK/2,
                ifelse(v > 0,
                    2*(v/PQBRAK)^2,
                    0),
                ifelse(v < PQBRAK,
                    1 - 2*((v - PQBRAK)/PQBRAK)^2,
                    1)),
            # Characteristic 2: Trigonometric behavior below PQBRAK
            ifelse(v < PQBRAK,
                a0 + a1*cos(v*wp) + b1*sin(v*wp),
                1))

        kI ~ ifelse(characteristic == 1,
            1,  # Always 1 for characteristic 1
            # Characteristic 2: Exponential behavior below 0.5 pu
            ifelse(v < 0.5,
                a2*b2*v^(b2 - 1)*exp(-a2*v^b2),
                1))

        # Main ZIP load equations with load variation (from OpenIPSL Load model)
        # kP can get to Inf for S_P_re + d_P = 0 so we limit it to Inf-ε
        kI*S_I_re*v + S_Y_re*v^2 + min(prevfloat(Inf), kP)*(S_P_re + d_P) ~ pvr*pir + pvi*pii
        kI*S_I_im*v + S_Y_im*v^2 + min(prevfloat(Inf), kP)*(S_P_im + d_Q) ~ (-pvr*pii) + pvi*pir
    end
end
