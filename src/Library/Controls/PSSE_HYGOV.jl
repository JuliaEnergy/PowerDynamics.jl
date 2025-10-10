# PSSE Hydro Turbine Governor Components
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.

@mtkmodel PSSE_HYGOV begin
    @parameters begin
        # Governor parameters
        R=0.05, [description="Permanent droop gain [pu]"]
        r=0.3, [description="Temporary droop gain [pu]"]
        T_r=5, [description="Governor time constant [s]"]
        T_f=0.05, [description="Filter time constant [s]"]
        T_g=0.5, [description="Servo time constant [s]"]
        VELM=0.2, [description="Gate open/close velocity limit [pu/s]"]
        G_MAX=0.9, [description="Maximum gate limit [pu]"]
        G_MIN=0, [description="Minimum gate limit [pu]"]

        # Hydraulic turbine parameters
        T_w=1.25, [description="Water time constant [s]"]
        A_t=1.2, [description="Turbine gain [pu]"]
        D_turb=0.2, [description="Turbine damping [pu]"]
        q_NL=0.08, [description="Water flow at no load [pu]"]
        h0=1, [guess=0.1, description="Water head initial value [pu]", bounds=(0, Inf)]

        # Initialization parameters (determined during initialization)
        n_ref, [guess=1, description="Speed reference [pu]", bounds=(0,Inf)]
    end

    # Protected parameters calculated from initial conditions (like OpenIPSL)
    # These follow the OpenIPSL initialization logic:
    # P_m0 := PMECH0
    # q0 := P_m0/(A_t*h0) + q_NL
    # g0 := q0/sqrt(h0)
    # c0 := g0
    # nref := R*c0

    @components begin
        # Active inputs/outputs
        SPEED_in = RealInput()   # Machine speed deviation from nominal [pu]
        PMECH_out = RealOutput() # Turbine mechanical power [pu]

        # Building block components
        filter = SimpleLag(K=1, T=T_f, default=0)
        temp_droop = SimpleLead(K=r*T_r, T=T_r, guess=0)
        position_limiter = LimIntegrator(K=1, outMin=G_MIN, outMax=G_MAX, guess=0.1)
        servo = SimpleLag(K=1, T=T_g, guess=0.1)
        water_integrator = LimIntegrator(K=1/T_w, outMin=-Inf, outMax=Inf, guess=0.1)
    end

    @variables begin
        # Governor system variables
        speed_input(t), [description="Speed input signal [pu]"]
        droop_feedback(t), [description="Permanent droop feedback [pu]"]
        governor_error(t), [description="Governor error signal [pu]"]
        velocity_limited_signal(t), [description="Velocity limited signal [pu/s]"]
        desired_gate_position(t), [description="Desired gate position [pu]"]
        gate_position(t), [description="Actual gate position [pu]"]

        # Hydraulic system variables
        flow_gate_ratio(t), [description="Flow to gate ratio"]
        flow_gate_ratio_squared(t), [description="(Q/G)^2"]
        water_flow(t), [description="Water flow Q [pu]"]
        flow_minus_noload(t), [description="Flow minus no-load Q-qNL [pu]"]
        turbine_head(t), [description="Turbine head H [pu]"]
        head_flow_product(t), [description="H*(Q-qNL) [pu]"]
        base_power(t), [description="A_t * H * (Q-qNL) [pu]"]
        damping_term(t), [description="D_turb * G * SPEED [pu]"]
    end

    @equations begin
        # Governor system (following OpenIPSL signal flow)
        speed_input ~ SPEED_in.u

        # Permanent droop feedback (R * gate_position)
        droop_feedback ~ R * desired_gate_position

        # Main error calculation (add: n_ref - (SPEED + R*c))
        governor_error ~ n_ref - (speed_input + droop_feedback)

        # Filter (SimpleLag1)
        filter.in ~ governor_error

        # Temporary droop compensation (simpleLead)
        temp_droop.in ~ filter.out

        # Velocity limiting
        velocity_limited_signal ~ clamp(temp_droop.out, -VELM, VELM)

        # Position limiter (LimIntegrator)
        position_limiter.in ~ velocity_limited_signal
        desired_gate_position ~ position_limiter.out

        # Servo motor
        servo.in ~ desired_gate_position
        gate_position ~ servo.out

        # Hydraulic turbine system (following OpenIPSL connections)
        # Division: Q/G
        flow_gate_ratio ~ water_flow / gate_position

        # Product: (Q/G)^2
        turbine_head ~ flow_gate_ratio^2

        # Water integrator: dQ/dt = (h0 - (Q/G)^2)/T_w
        water_integrator.in ~ h0 - turbine_head
        water_flow ~ water_integrator.out

        # Add3: Q - qNL
        # flow_minus_noload ~ water_flow - q_NL

        # Gain6: A_t * H * (Q - qNL)
        base_power ~ A_t * turbine_head * (water_flow - q_NL)

        # Damping term: D_turb * G * SPEED
        damping_term ~ D_turb * gate_position * speed_input

        # Final power: base_power - damping_term
        PMECH_out.u ~ base_power - damping_term
    end
end
