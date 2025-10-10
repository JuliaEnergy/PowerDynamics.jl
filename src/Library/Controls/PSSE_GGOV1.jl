# PSSE Turbine Governor Components
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.

@mtkmodel _GGOV1_AccelerationLimiter begin
    @structural_parameters begin
        Ka # Acceleration limiter gain [pu]
        Ta # Acceleration limiter time constant [s]
        DELT # Time step used in simulation [s]
    end

    @components begin
        speed_derivative = Derivative(K=1, T=Ta)
    end

    @variables begin
        SPEED(t), [input=true, description="Machine speed deviation from nominal [pu]"]
        ASET(t), [input=true, description="Acceleration limiter setpoint [pu/s]"]
        FSR(t), [input=true, description="Governor output signal [pu]"]
        FSRA(t), [output=true, description="Acceleration controller output [pu]"]
        # Internal signals
        error_signal(t), [description="ASET minus speed derivative [pu/s]"]
        acceleration_output(t), [description="Scaled acceleration output [pu]"]
    end
    @equations begin
        speed_derivative.in ~ SPEED
        FSRA ~ FSR + Ka * DELT * (ASET - speed_derivative.out)
    end
end

@mtkmodel _GGOV1_LoadLimiter begin
    @structural_parameters begin
        Kturb # Turbine gain
        Kpload # Load limiter proportional gain for PI controller
        Kiload # Load limiter integral gain for PI controller
        Dm # Mechanical damping coefficient
        Wfnl # No load fuel flow
    end

    @variables begin
        LDREF(t), [input=true, description="Load limiter reference value [pu]"]
        TEXM(t), [input=true, description="Measured exhaust temperature [pu]"]
        PELEC(t), [input=true, description="Machine electrical power [pu]"]
        FSRT(t), [output=true, description="Controller output [pu]"]
        # Internal state and signals
        integral_state(t), [description="Integrator state"]
        temp_error(t), [description="Temperature error signal"]
    end

    @equations begin
        # Temperature error calculation
        temp_error ~ (LDREF / Kturb + Wfnl) - TEXM

        # PI controller with direct integration
        Dt(integral_state) ~ Kiload * temp_error

        # Output with PI terms and upper limiting
        FSRT ~ min(integral_state + Kpload * temp_error, 1.0)
    end
end

@mtkmodel _GGOV1_PIDGovernor begin
    @structural_parameters begin
        Rselect # Feedback signal selector (1/-1/-2/0)
        R # Permanent droop
        T_pelec # Electrical power transducer time constant
        maxerr # Maximum speed error
        minerr # Minimum speed error
        Kpgov # Governor proportional gain
        Kigov # Governor integral gain
        Kdgov # Governor derivative gain
        Tdgov # Governor derivative time constant
        Kturb # Turbine gain (for interface compatibility)
        Kimw # Power controller gain
        db # Speed governor deadband
        Wfnl # No load fuel flow (for interface compatibility)
        Dm # Mechanical damping coefficient (for interface compatibility)
    end

    @components begin
        power_transducer = SimpleLag(K=1, T=T_pelec, guess=nothing)
        power_controller = LimIntegrator(K=Kimw, outMin=-1.1*R, outMax=1.1*R, guess=nothing)
        speed_derivative = Derivative(K=Kdgov, T=Tdgov)
        deadband = DeadZone(uMax=db, uMin=-db)
    end

    @variables begin
        # Inputs
        PELEC(t), [input=true, description="Machine electrical power [pu]"]
        PMW_SET(t), [input=true, description="Supervisory load controller setpoint [pu]"]
        P_REF(t), [input=true, description="Power reference [pu]"]
        SPEED(t), [input=true, description="Machine speed deviation [pu]"]
        VSTROKE(t), [input=true, description="Valve stroke [pu]"]
        GOVOUT1(t), [input=true, description="Governor output before limiter [pu]"]
        FSRN(t), [output=true, description="Governor output [pu]"]

        # Internal states and signals
        pid_integral_state(t), [description="PID integral state"]

        # Internal signals
        selected_feedback(t), [description="R_select output"]
        speed_error(t), [description="Speed error before deadband"]
        limited_error(t), [description="Limited speed error"]
    end

    @equations begin
        # Power transducer
        power_transducer.in ~ PELEC

        # Power controller (LimIntegrator with limits ±1.1*R)
        power_controller.in ~ PMW_SET - power_transducer.out

        # R_select logic (feedback signal selection) - structural parameter allows if/else
        if Rselect == 1
            selected_feedback ~ power_transducer.out
        elseif Rselect == -1
            selected_feedback ~ VSTROKE
        elseif Rselect == -2
            selected_feedback ~ GOVOUT1
        else
            selected_feedback ~ 0.0
        end

        # Speed error calculation
        speed_error ~ -SPEED + P_REF + power_controller.out - R * selected_feedback

        # Deadband
        deadband.in ~ speed_error

        # Error limiting
        limited_error ~ clamp(deadband.out, minerr, maxerr)

        # PID controller
        speed_derivative.in ~ limited_error
        Dt(pid_integral_state) ~ Kigov * limited_error

        # PID output
        FSRN ~ Kpgov * limited_error + pid_integral_state + speed_derivative.out
    end
end

@mtkmodel _GGOV1_Turbine begin
    @structural_parameters begin
        Flag # Switch for fuel source (0/1)
        Tact # Actuator time constant
        Kturb # Turbine gain
        Tb # Turbine lag time constant
        Tc # Turbine lead time constant
        Teng # Transport lag (must be 0)
        Tfload # Load limiter time constant
        Dm # Mechanical damping coefficient
        Vmax # Maximum valve position
        Vmin # Minimum valve position
        Ropen # Maximum valve opening rate
        Rclose # Maximum valve closing rate (negative)
        Tsa # Temperature detection lead time
        Tsb # Temperature detection lag time
        DELT # Time step (for interface compatibility)
        Wfnl # No load fuel flow
    end

    begin
        # Validation check for unsupported delay
        Teng > 0 && error("Teng > 0 not supported!")
    end

    @components begin
        turbine_dynamics = LeadLag(K=1, T1=Tc, T2=Tb, guess=nothing)
        temp_leadlag = LeadLag(K=1, T1=Tsa, T2=Tsb, guess=nothing)
        temp_lag = SimpleLag(K=1, T=Tfload, guess=nothing)
    end

    @variables begin
        # Inputs
        FSR(t), [input=true, description="Governor output [pu]"]
        SPEED(t), [input=true, description="Machine speed deviation [pu]"]
        PELEC(t), [input=true, description="Machine electrical power [pu]"]

        # Outputs
        PMECH(t), [output=true, description="Turbine mechanical power [pu]"]
        VSTROKE(t), [output=true, description="Valve position [pu]"]
        TEXM(t), [output=true, description="Measured exhaust temperature [pu]"]

        # Internal states
        valve_integrator(t), [description="Valve actuator integrator state"]

        # Internal signals
        speed_offset(t), [description="Speed offset by 1.0"]
        flag_output(t), [description="Flag logic output"]
        dm_output(t), [description="Dm_select logic output"]
        valve_error(t), [description="FSR - valve_position"]
        valve_rate(t), [description="Rate-limited valve error"]
        fuel_flow(t), [description="Flag_output * valve_position"]
        turbine_input(t), [description="Kturb * (fuel_flow - Wfnl)"]
        temp_input(t), [description="Temperature path input"]
    end

    @equations begin
        # Flag logic (structural allows if/else)
        if Flag == 1
            flag_output ~ SPEED + 1
        else
            flag_output ~ 1.0
        end

        # Dm_select logic
        dm_output ~ ifelse(Dm ≥ 0, (SPEED + 1)*Dm, (SPEED + 1)^Dm)

        # Actuator dynamics (rate and position limited)
        valve_error ~ FSR - VSTROKE
        valve_rate ~ clamp(valve_error / Tact, Rclose, Ropen)
        Dt(valve_integrator) ~ valve_rate
        VSTROKE ~ clamp(valve_integrator, Vmin, Vmax)

        # Fuel flow calculation
        fuel_flow ~ flag_output * VSTROKE
        turbine_input ~ Kturb * (fuel_flow - Wfnl)

        # Turbine dynamics (no delay)
        turbine_dynamics.in ~ turbine_input

        # Temperature path
        # XXX: DEVEATION from OpenIpsl: only multipy when Dm < 0
        temp_input ~ fuel_flow * ifelse(Dm ≥ 0, 1, dm_output)
        temp_leadlag.in ~ temp_input
        temp_lag.in ~ temp_leadlag.out

        # Outputs
        # XXX: DEVIATIOn from OpenIPSL: only substract when Dm ≥ 0
        PMECH ~ turbine_dynamics.out - ifelse(Dm ≥ 0, dm_output, 0)
        TEXM ~ temp_lag.out
    end
end

@mtkmodel PSSE_GGOV1_EXPERIMENTAL begin
    @structural_parameters begin
        Flag=1 # Switch for fuel source characteristic (0/1)
        Rselect=1 # Feedback signal for governor droop (1/-1/-2/0)
        Teng=0 # Transport lag time constant for diesel engine [s]
    end
    @parameters begin
        # Governor parameters
        # Rselect as structural parameter
        _Rselect_static=Rselect, [description="Static copy of Rselect for internal use"]
        R=0.04, [description="Permanent droop [pu]"]
        T_pelec=1, [description="Electrical power transducer time constant [s]"]
        maxerr=0.05, [description="Maximum value for speed error signal [pu]"]
        minerr=-0.05, [description="Minimum value for speed error signal [pu]"]
        Kpgov=10, [description="Governor proportional gain [pu]"]
        Kigov=2, [description="Governor integral gain [pu]"]
        Kdgov=0, [description="Governor derivative gain [pu]"]
        Tdgov=1, [description="Governor derivative controller time constant [s]"]
        Kimw=0, [description="Power controller (reset) gain [pu]"]
        db=0, [description="Speed governor deadband [pu]"]

        # Turbine parameters
        # FLAG as structural parameter
        _Flag_static = Flag, [description="Static copy of Flag for internal use"]
        Tact=0.5, [description="Actuator time constant [s]"]
        Kturb=1.5, [description="Turbine gain [pu]"]
        Tb=0.1, [description="Turbine lag time constant [s]"]
        Tc=0, [description="Turbine lead time constant [s]"]
        # Teng as structural parameter to throw error
        # Teng=0, [description="Transport lag time constant for diesel engine [s]"]
        Tsa=4, [description="Temperature detection lead time constant [s]"]
        Tsb=5, [description="Temperature detection lag time constant [s]"]

        # Load limiter parameters
        Tfload=3, [description="Load Limiter time constant [s]"]
        Kpload=2, [description="Load limiter proportional gain for PI controller [pu]"]
        Kiload=0.67, [description="Load limiter integral gain for PI controller [pu]"]
        Ldref, [description="Load limiter reference value [pu]"]

        # Acceleration limiter parameters
        Ka=10, [description="Acceleration limiter gain [pu]"]
        Ta=0.1, [description="Acceleration limiter time constant [s]"]
        Aset=0.1, [description="Acceleration limiter setpoint [pu/s]"]

        # Valve and actuator limits
        Vmax=1, [description="Maximum valve position limit [pu]"]
        Vmin=0.15, [description="Minimum valve position limit [pu]"]
        Ropen=0.1, [description="Maximum valve opening rate [pu/s]"]
        Rclose=-0.1, [description="Maximum valve closing rate [pu/s]"]

        # Common parameters
        Dm=0, [description="Mechanical damping coefficient [pu]"]
        Wfnl=0.2, [description="No load fuel flow [pu]"]
        DELT=0.005, [description="Time step used in simulation [s]"]

        # Initialization parameters (determined during initialization)
        Pref, [description="Power reference setpoint [pu]"]
        Pmwset, [description="Supervisory load controller setpoint [pu]"]
    end

    @components begin
        pid_governor = _GGOV1_PIDGovernor(
            Rselect=Rselect, R=R, T_pelec=T_pelec, maxerr=maxerr, minerr=minerr,
            Kpgov=Kpgov, Kigov=Kigov, Kdgov=Kdgov, Tdgov=Tdgov, Kturb=Kturb,
            Kimw=Kimw, db=db, Wfnl=Wfnl, Dm=Dm
        )
        load_limiter = _GGOV1_LoadLimiter(
            Kturb=Kturb, Kpload=Kpload, Kiload=Kiload, Dm=Dm, Wfnl=Wfnl
        )
        accel_limiter = _GGOV1_AccelerationLimiter(
            Ka=Ka, Ta=Ta, DELT=DELT
        )
        turbine = _GGOV1_Turbine(
            Flag=Flag, Tact=Tact, Kturb=Kturb, Tb=Tb, Tc=Tc, Teng=Teng,
            Tfload=Tfload, Dm=Dm, Vmax=Vmax, Vmin=Vmin, Ropen=Ropen,
            Rclose=Rclose, Tsa=Tsa, Tsb=Tsb, DELT=DELT, Wfnl=Wfnl
        )

        SPEED_in = RealInput() # Input: Machine speed deviation from nominal [pu]
        PELEC_in = RealInput() # Input: Machine electrical power [pu]
        PMECH_out = RealOutput() # Output: Turbine mechanical power [pu]
    end

    @variables begin
        # Internal signals
        min_select_out(t), [description="Minimum of three controller outputs [pu]"]
        fsr_limited(t), [description="Position-limited governor output [pu]"]
    end

    @equations begin
        # Input distribution to all submodels
        pid_governor.SPEED ~ SPEED_in.u
        pid_governor.PELEC ~ PELEC_in.u
        pid_governor.P_REF ~ Pref
        pid_governor.PMW_SET ~ Pmwset

        load_limiter.PELEC ~ PELEC_in.u
        load_limiter.LDREF ~ Ldref

        accel_limiter.SPEED ~ SPEED_in.u
        accel_limiter.ASET ~ Aset

        turbine.SPEED ~ SPEED_in.u
        turbine.PELEC ~ PELEC_in.u
        turbine.FSR ~ fsr_limited

        # Min select logic (direct implementation)
        min_select_out ~ min(pid_governor.FSRN, load_limiter.FSRT, accel_limiter.FSRA)

        # Position limiter for FSR
        fsr_limited ~ clamp(min_select_out, Vmin, Vmax)

        # FSR feedback distribution
        accel_limiter.FSR ~ fsr_limited
        pid_governor.GOVOUT1 ~ fsr_limited
        pid_governor.VSTROKE ~ turbine.VSTROKE

        # Temperature feedback loop
        load_limiter.TEXM ~ turbine.TEXM

        # Final output
        PMECH_out.u ~ turbine.PMECH
    end
end
