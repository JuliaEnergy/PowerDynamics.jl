# PSSE Excitation System Components
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.

@mtkmodel RotatingExciter begin
    @structural_parameters begin
        T_E # Exciter time constant
        K_E # Exciter field gain
        E_1 # Exciter saturation point 1
        E_2 # Exciter saturation point 2
        S_EE_1 # Saturation at E_1
        S_EE_2 # Saturation at E_2
    end
    @variables begin
        I_C(t), [description="Control input current", input=true]
        EFD(t), [guess=1, description="Exciter field voltage output", output=true]
        # intermediate variables
        feedback(t)
    end
    @equations begin
        Dt(EFD) ~ 1/T_E * (I_C - feedback)
        feedback ~ K_E * EFD + (EFD * PSSE_QUAD_SE(EFD, S_EE_1, S_EE_2, E_1, E_2))
    end
end

@mtkmodel PSSE_IEEET1 begin
    @parameters begin
        # Public parameters with OpenIPSL defaults
        T_R=0.02, [description="Regulator input filter time constant [s]"]
        K_A=200, [description="Regulator output gain [pu]"]
        T_A=0.001, [description="Regulator output time constant [s]"]
        V_RMAX=2, [description="Maximum regulator output [pu]"]
        V_RMIN=-2, [description="Minimum regulator output [pu]"]
        K_E=0.1, [description="Exciter field proportional constant [pu]"]
        T_E=0.55, [description="Exciter field time constant [s]"]
        K_F=0.06, [description="Rate feedback excitation system stabilizer gain [pu]"]
        T_F=1, [description="Rate feedback time constant [s]"]
        E_1=2.85, [description="Exciter output voltage for saturation factor S_E(E_1) [pu]"]
        S_EE_1=0.3, [description="Exciter saturation factor at exciter output voltage E1 [pu]"]
        E_2=3.8, [description="Exciter output voltage for saturation factor S_E(E_2) [pu]"]
        S_EE_2=0.6, [description="Exciter saturation factor at exciter output voltage E2 [pu]"]

        # Free initialization parameter
        V_REF, [guess=1, description="Voltage reference setpoint [pu]"]
    end

    @components begin
        # Active inputs/outputs
        ECOMP_in = RealInput()    # Terminal voltage measurement input
        # EFD0_in = RealInput()     # Initial field voltage for initialization
        EFD_out = RealOutput()    # Field voltage output to generator

        # Unused inputs (commented out but kept for reference)
        # VOTHSG_in = RealInput()   # Other signal input (typically zero)
        # VUEL_in = RealInput()     # Under-excitation limiter input
        # VOEL_in = RealInput()     # Over-excitation limiter input
        # XADIFD_in = RealInput()   # Machine field current input

        # Building block components
        transducer = SimpleLag(K=1, T=T_R)
        amplifier = SimpleLagLim(K=K_A, T=T_A, outMin=V_RMIN, outMax=V_RMAX)
        derivative_lag = Derivative(K=K_F, T=T_F)
        exciter = RotatingExciter(T_E=T_E, K_E=K_E, E_1=E_1, E_2=E_2, S_EE_1=S_EE_1, S_EE_2=S_EE_2)
    end

    @variables begin
        # Internal signal variables
        voltage_error(t), [description="Voltage error signal [pu]"]
        sum_signal(t), [description="Summed signal before amplifier [pu]"]
        derivative_feedback(t), [description="Derivative feedback signal [pu]"]
    end

    @equations begin
        # Input processing
        transducer.in ~ ECOMP_in.u

        # Voltage error calculation
        voltage_error ~ V_REF - transducer.out

        # Signal summing (3-input addition: voltage_error + VOTHSG - derivative_feedback)
        # Note: VOTHSG and limiter signals are omitted (typically zero)
        sum_signal ~ voltage_error - derivative_feedback

        # Amplifier with limiting
        amplifier.in ~ sum_signal

        # Exciter
        exciter.I_C ~ amplifier.out

        # Derivative feedback
        derivative_lag.in ~ exciter.EFD
        derivative_feedback ~ derivative_lag.out

        # Output connection
        EFD_out.u ~ exciter.EFD
    end
end

@mtkmodel PSSE_EXST1 begin
    @structural_parameters begin
        vothsg_input = false  # Other signal input (rarely used)
        vuel_input = false    # Under-excitation limiter input (rarely used)
        voel_input = false    # Over-excitation limiter input (rarely used)
    end

    @parameters begin
        # Fixed parameters with OpenIPSL defaults
        T_R=0.02, [description="Regulator input filter time constant [s]"]
        V_IMAX=0.2, [description="Maximum voltage error (regulator input) [pu]"]
        V_IMIN=0, [description="Minimum voltage error (regulator input) [pu]"]
        T_C=1, [description="Regulator numerator (lead) time constant [s]"]
        T_B=1, [description="Regulator denominator (lag) time constant [s]"]
        K_A=80, [description="Voltage regulator gain [pu]"]
        T_A=0.05, [description="Voltage regulator time constant [s]"]
        V_RMAX=8, [description="Maximum exciter output [pu]"]
        V_RMIN=-3, [description="Minimum exciter output [pu]"]
        K_C=0.2, [description="Rectifier loading factor proportional to commutating reactance [pu]"]
        K_F=0.1, [description="Rate feedback gain [pu]"]
        T_F=1, [description="Rate feedback time constant [s]"]

        # Optional input default values (when inputs are disabled)
        if !vothsg_input
            VOTHSG_set=0.0, [description="Other signal setpoint [pu]"]
        end
        if !vuel_input
            VUEL_set=0.0, [description="Under-excitation limiter setpoint [pu]"]
        end
        if !voel_input
            VOEL_set=0.0, [description="Over-excitation limiter setpoint [pu]"]
        end

        # Free initialization parameter
        V_REF, [guess=1, description="Voltage reference setpoint [pu]"]
    end

    @components begin
        # Always required input/output interfaces
        ECOMP_in = RealInput()     # Terminal voltage measurement
        EFD_out = RealOutput()     # Field voltage output
        XADIFD_in = RealInput()    # Machine field current

        # Optional auxiliary inputs (conditional on structural parameters)
        if vothsg_input
            VOTHSG_in = RealInput()    # Other signal input
        end
        if vuel_input
            VUEL_in = RealInput()      # Under-excitation limiter
        end
        if voel_input
            VOEL_in = RealInput()      # Over-excitation limiter
        end

        # Building block components
        transducer = SimpleLag(K=1, T=T_R, guess=1)
        leadlag = LeadLag(K=1, T1=T_C, T2=T_B)
        amplifier = SimpleLag(K=1, T=T_A)
        derivative_feedback = Derivative(K=K_F, T=T_F)
    end

    @variables begin
        # Signal processing variables
        voltage_error(t), [description="Voltage error signal [pu]"]
        sum_signal(t), [description="Summed signal (inlcuding feedback) [pu]"]
        voltage_limited(t), [description="Limited lead lag input [pu]"]
        EFD_unlimited(t), [description="Amplifier output before limiting [pu]"]
        vr_max_limit(t), [description="Upper voltage limit with rectifier drop [pu]"]
        vr_min_limit(t), [description="Lower voltage limit with rectifier drop [pu]"]
    end

    @equations begin
        # Input transducer
        transducer.in ~ ECOMP_in.u

        # Voltage error calculation
        voltage_error ~ V_REF - transducer.out

        # Signal summing: lead-lag output + optional inputs - derivative feedback
        sum_signal ~ voltage_error +
                     (vothsg_input ? VOTHSG_in.u : VOTHSG_set) +
                     (vuel_input ? VUEL_in.u : VUEL_set) +
                     (voel_input ? VOEL_in.u : VOEL_set) -
                     derivative_feedback.out

        voltage_limited ~ clamp(sum_signal, V_IMIN, V_IMAX)

        # Lead-lag compensator
        leadlag.in ~ voltage_limited

        # Amplifier with limiting
        amplifier.in ~ K_A * leadlag.out
        EFD_unlimited ~ amplifier.out

        # Derivative feedback
        derivative_feedback.in ~ EFD_unlimited

        # Rectifier commutation voltage drop limits
        vr_max_limit ~ ECOMP_in.u * V_RMAX - K_C * XADIFD_in.u
        vr_min_limit ~ ECOMP_in.u * V_RMIN - K_C * XADIFD_in.u

        # Final output with rectifier commutation limiting
        EFD_out.u ~ clamp(EFD_unlimited, vr_min_limit, vr_max_limit)
    end
end

@mtkmodel PSSE_ESST4B begin
    @structural_parameters begin
        vothsg_input = false  # Other signal input (rarely used)
        vuel_input = false    # Under-excitation limiter input (rarely used)
        voel_input = false    # Over-excitation limiter input (rarely used)
    end

    @parameters begin
        # Fixed parameters with OpenIPSL defaults
        T_R=0.3, [description="Regulator input filter time constant [s]"]
        K_PR=2.97, [description="Voltage regulator proportional gain [pu]"]
        K_IR=2.97, [description="Voltage regulator integral gain [pu]"]
        V_RMAX=1, [description="Maximum regulator output [pu]"]
        V_RMIN=-0.87, [description="Minimum regulator output [pu]"]
        T_A=0.01, [description="Thyristor bridge time constant [s]"]
        K_PM=1, [description="Current regulator proportional gain [pu]"]
        K_IM=0.2, [description="Current regulator integral gain [pu]"]
        V_MMAX=1, [description="Maximum current regulator output [pu]"]
        V_MMIN=-0.87, [description="Minimum current regulator output [pu]"]
        K_G=0.1, [description="Feedback gain [pu]"]
        K_P=6.73, [description="Potential circuit gain [pu]"]
        K_I=0.1, [description="Current circuit gain [pu]"]
        V_BMAX=8.41, [description="Maximum exciter voltage [pu]"]
        K_C=0.1, [description="Rectifier loading factor [pu]"]
        X_L=0, [description="Reactance associated with potential source [pu]"]
        THETAP=0, [description="Potential circuit phase angle [rad]"]

        # Optional input default values (when inputs are disabled)
        if !vothsg_input
            VOTHSG_set=0.0, [description="Other signal setpoint [pu]"]
        end
        if !vuel_input
            VUEL_set=0.0, [description="Under-excitation limiter setpoint [pu]"]
        end
        if !voel_input
            VOEL_set=1e10, [description="Over-excitation limiter setpoint [pu]"]
        end

        # Free initialization parameter
        V_REF, [guess=1, description="Voltage reference setpoint [pu]"]
    end

    @components begin
        # Always required input/output interfaces
        ECOMP_in = RealInput()     # Terminal voltage measurement
        EFD_out = RealOutput()     # Field voltage output
        XADIFD_in = RealInput()    # Machine field current

        # Terminal inputs (from PSSE_BaseMachine outputs)
        TERM_VR_in = RealInput()   # Terminal voltage real part
        TERM_VI_in = RealInput()   # Terminal voltage imag part
        TERM_IR_in = RealInput()   # Terminal current real part
        TERM_II_in = RealInput()   # Terminal current imag part

        # Optional auxiliary inputs (conditional on structural parameters)
        if vothsg_input
            VOTHSG_in = RealInput()    # Other signal input
        end
        if vuel_input
            VUEL_in = RealInput()      # Under-excitation limiter
        end
        if voel_input
            VOEL_in = RealInput()      # Over-excitation limiter
        end

        # Building block components
        transducer = SimpleLag(K=1, T=T_R)
        voltage_int = LimIntegrator(K=K_IR, outMin=V_RMIN/K_PR, outMax=V_RMAX/K_PR)
        thyristor = SimpleLag(K=1, T=T_A)
        current_int = LimIntegrator(K=K_IM, outMin=V_MMIN/K_PM, outMax=V_MMAX/K_PM)
        rectifier = RectifierCommutationVoltageDrop(K_C=K_C)
    end

    @variables begin
        # Signal processing variables
        voltage_error(t), [description="Voltage error signal [pu]"]
        vr_sum(t), [description="Voltage regulator input sum [pu]"]
        vr_prop(t), [description="Voltage regulator proportional output [pu]"]
        vr_out(t), [guess=1, description="Voltage regulator output [pu]"]
        va_out(t), [description="Thyristor bridge output [pu]"]
        current_error(t), [description="Current error with feedback [pu]"]
        vm_prop(t), [guess=0, description="Current regulator proportional output [pu]"]
        vm_out(t), [guess=1, description="Current regulator output [pu]"]
        vm_limited(t), [description="LV_GATE limited output [pu]"]
        vb_signal(t), [description="Exciter output (limited) [pu]"]
        VE(t), [description="adapted terminal voltage magnitude for rectifier"]
    end
    begin
        V_T = TERM_VR_in.u + im*TERM_VI_in.u
        I_T = TERM_IR_in.u + im*TERM_II_in.u
        K_P_comp = K_P*cos(THETAP) + im*K_P*sin(THETAP)
        ve_term = simplify(abs(K_P_comp*V_T + im*(K_I + K_P_comp*X_L)*I_T))
    end

    @equations begin
        # Input transducer
        transducer.in ~ ECOMP_in.u

        # Voltage error calculation
        voltage_error ~ V_REF - transducer.out

        # 3-input sum with conditional switching
        vr_sum ~ voltage_error +
                  (vothsg_input ? VOTHSG_in.u : VOTHSG_set) +
                  (vuel_input ? VUEL_in.u : VUEL_set)

        # Voltage PI controller
        voltage_int.in ~ vr_sum
        vr_prop ~ K_PR * vr_sum
        vr_out ~ clamp(vr_prop + voltage_int.out, V_RMIN, V_RMAX)

        # Thyristor bridge
        thyristor.in ~ vr_out
        va_out ~ thyristor.out

        # Current error with feedback
        current_error ~ va_out - K_G * EFD_out.u

        # Current PI controller
        current_int.in ~ current_error
        vm_prop ~ K_PM * current_error
        vm_out ~ clamp(vm_prop + current_int.out, V_MMIN, V_MMAX)

        # LV_GATE with conditional VOEL
        vm_limited ~ min(vm_out, voel_input ? VOEL_in.u : VOEL_set)

        # Rectifier commutation voltage drop
        VE ~ ve_term
        rectifier.V_EX ~ VE
        rectifier.XADIFD ~ XADIFD_in.u

        # Final limiting with V_BMAX
        vb_signal ~ clamp(rectifier.EFD, -Inf, V_BMAX)

        # Output connection
        EFD_out.u ~ vb_signal * vm_limited
    end
end
