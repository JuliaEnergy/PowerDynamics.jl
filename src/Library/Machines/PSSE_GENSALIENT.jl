# PSSE GENSALIENT Model - Unified Salient Pole Generator (GENSAL/GENSAE)
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.
#
# Description: Unified salient pole generator model with configurable saturation

#=
On structural parameters:
SE is a structural parameter allowing different saturation functions to be passed in

Key differences from GENROUND (round rotor):
- Only 3 state variables: Epq, PSIkd, PSIppq (no Epd, PSIkq)
- No q-axis transient dynamics (no Tpq0 parameter)
- Simpler equations for salient pole behavior
- Different field current calculations
=#

@mtkmodel PSSE_GENSALIENT begin
    @structural_parameters begin
        pmech_input = true
        efd_input = true
        SE  # Saturation function - can be PSSE_QUAD_SE or PSSE_EXP_SE
    end
    @extend PSSE_BaseMachine()

    @components begin
        if pmech_input
            PMECH_in = RealInput() # Turbine mechanical power
        end
        if efd_input
            EFD_in = RealInput()   # Generator main field voltage
        end
    end

    @parameters begin
        # we cannot @extend a model with structural parameters, so we need to define it in concrete mode
        if !pmech_input
            pmech_set, [guess=1, description="mechanical power setpoint [pu]"]
        end
        if !efd_input
            efd_set, [guess=1, description="field voltage setpoint [pu]"]
        end
        # Note: No Xpq and Tpq0 parameters for salient pole models
        # q-axis transient reactance equals q-axis reactance for salient pole
    end

    @variables begin
        # State variables (3 states for salient pole, vs 4 for round rotor)
        Epq(t), [guess=1, description="q-axis voltage behind transient reactance [pu]"]
        PSIkd(t), [guess=1, description="d-axis rotor flux linkage [pu]"]
        PSIppq(t), [guess=0, description="q-axis subtransient flux linkage [pu]"]

        # Algebraic variables
        PSId(t), [guess=1, description="d-axis flux linkage [pu]"]
        PSIq(t), [guess=0, description="q-axis flux linkage [pu]"]
        PSIppd(t), [guess=1, description="d-axis subtransient flux linkage [pu]"]
        XadIfd(t), [guess=1, description="d-axis machine field current [pu]"]
    end

    begin
        # Constants (from OpenIPSL GENSAL/GENSAE lines 77-80, 86-89)
        K1d = (Xpd - Xppd)*(Xd - Xpd)/(Xpd - Xl)^2
        K2d = (Xpd - Xl)*(Xppd - Xl)/(Xpd - Xppd)
        K3d = (Xppd - Xl)/(Xpd - Xl)
        K4d = (Xpd - Xppd)/(Xpd - Xl)
    end

    @equations begin
        # Input switching (we cannot extend a model with structural parameters)
        pmech ~ pmech_input ? PMECH_in.u : pmech_set
        efd ~ efd_input ? EFD_in.u : efd_set

        # State equations (from OpenIPSL GENSAL lines 91-93, GENSAE lines 100-102)
        # Note: Only 3 differential equations for salient pole models
        Dt(Epq) ~ 1/Tpd0*(efd - XadIfd)
        Dt(PSIkd) ~ 1/Tppd0*(Epq - PSIkd - (Xpd - Xl)*id)

        # Q-axis subtransient equation differs between GENSAL and GENSAE
        # GENSAL (quadratic): simpler equation without saturation
        # GENSAE (exponential): includes saturation term
        # Use ifelse to switch based on saturation function type
        if SE==PSSE_QUAD_SE
            # GENSAL equation (line 93)
            Dt(PSIppq) ~ 1/Tppq0*((-PSIppq) + (Xq - Xppq)*iq)
        elseif SE==PSSE_EXP_SE
            # GENSAE equation (lines 102-107)
            Dt(PSIppq) ~ 1/Tppq0*((-PSIppq) + (Xq - Xppq)*iq - PSIppq*(Xq-Xl)/(Xd-Xl)*SE(sqrt(PSIppd*PSIppd + PSIppq*PSIppq), S10, S12, 1, 1.2))
        else
            error("Unsupported saturation function SE")
        end

        # Override electrical torque calculation (from OpenIPSL lines 103, 118)
        Te ~ PSId*iq - PSIq*id

        # Flux linkage equations (from OpenIPSL lines 94-96, 108-110)
        PSIppd ~ Epq*K3d + PSIkd*K4d
        PSId ~ PSIppd - Xppd*id
        PSIq ~ (-PSIppq) - Xppq*iq

        # Field current equations with configurable saturation
        # GENSAL (lines 97-102): uses saturation on Epq
        # GENSAE (lines 112-117): uses saturation on PSIpp magnitude
        if SE==PSSE_QUAD_SE
            # GENSAL field current equation
            XadIfd ~ K1d*(Epq - PSIkd - (Xpd - Xl)*id) + (Xd - Xpd)*id + (SE(Epq, S10, S12, 1, 1.2) + 1)*Epq
        elseif SE==PSSE_EXP_SE
            # GENSAE field current equation
            XadIfd ~ Epq + K1d*(Epq - PSIkd - (Xpd - Xl)*id) + (Xd - Xpd)*id + SE(sqrt(PSIppd*PSIppd + PSIppq*PSIppq), S10, S12, 1, 1.2)*PSIppd
        else
            error("Unsupported saturation function SE")
        end

        # Override voltage equations (from OpenIPSL lines 104-105, 119-120)
        ud ~ (-PSIq) - R_a*id
        uq ~ PSId - R_a*iq

        # Output connections (from OpenIPSL lines 87, 90, 96, 99)
        XADIFD_out.u ~ XadIfd
        ISORCE_out.u ~ XadIfd
    end
end

"""
    PSSE_GENSAL

This model is a port of the OpenIPSL [`Electrical.Machines.PSSE.GENSAL`](https://github.com/OpenIPSL/OpenIPSL/blob/fe8aa5c/OpenIPSL/Electrical/Machines/PSSE/GENSAL.mo) model,
maintaining the same mathematical formulation while adapting to PowerDynamics/ModelingToolkit conventions.

# Validation

Validated against the OpenIPSL SMIB testcase
[`Tests.Machines.PSSE.GENSAL`](https://github.com/OpenIPSL/OpenIPSL/blob/fe8aa5c/OpenIPSL/Tests/Machines/PSSE/GENSAL.mo).
See [validation plot](../assets/OpenIPSL_valid/GENSAL.png) generated by automatic validation script in `/test/OpenIPSL_test`.
"""
PSSE_GENSAL(; kwargs...) = PSSE_GENSALIENT(; SE=PSSE_QUAD_SE, kwargs...)

"""
    PSSE_GENSAE

This model is a port of the OpenIPSL [`Electrical.Machines.PSSE.GENSAE`](https://github.com/OpenIPSL/OpenIPSL/blob/fe8aa5c/OpenIPSL/Electrical/Machines/PSSE/GENSAE.mo) model,
maintaining the same mathematical formulation while adapting to PowerDynamics/ModelingToolkit conventions.

# Validation

Validated against the OpenIPSL SMIB testcase
[`Tests.Machines.PSSE.GENSAE`](https://github.com/OpenIPSL/OpenIPSL/blob/fe8aa5c/OpenIPSL/Tests/Machines/PSSE/GENSAE.mo).
See [validation plot](../assets/OpenIPSL_valid/GENSAE.png) generated by automatic validation script in `/test/OpenIPSL_test`.
"""
PSSE_GENSAE(; kwargs...) = PSSE_GENSALIENT(; SE=PSSE_EXP_SE, kwargs...)
