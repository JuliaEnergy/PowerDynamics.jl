# PSSE EXST1 Excitation System
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.

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
