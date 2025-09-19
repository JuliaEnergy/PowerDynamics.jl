# PSSE GENROU Model - Port of OpenIPSL.Electrical.Machines.PSSE.GENROU
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.
#
# Description: Round rotor generator model with quadratic saturation

"""
    PSSE_QUAD_SE(u, SE1, SE2, E1, E2)

Scaled Quadratic Saturation Function (PTI PSS/E).
Port of OpenIPSL.NonElectrical.Functions.PSSE_QUAD_SE
"""
function PSSE_QUAD_SE(u, SE1, SE2, E1, E2)
    if !(SE1 > 0.0 || SE1 < 0.0) || u <= 0.0
        return 0.0
    end

    # XXX: This is weird! The original code uses
    # parameter Real a=if not (SE2 > 0.0 or SE2 < 0.0) then sqrt(SE1*E1/(SE2*E2)) else 0;
    # which has the not operator so it is differnet from this julia function
    # maybe i missunderstand something about not or the or operator in modelica?
    a = if (SE2 > 0.0 || SE2 < 0.0)
        sqrt(SE1*E1/(SE2*E2))
    else
        0.0
    end

    A = E2 - (E1 - E2)/(a - 1)
    B = if abs(E1 - E2) < eps()
        0.0
    else
        SE2*E2*(a - 1)^2/(E1 - E2)^2
    end

    if u <= A
        return 0.0
    else
        return B*(u - A)^2/u
    end
end

#=
PSSE_QUAD_SSE uses if/else statements. We need to register it as a symbolic function
to block MTK from tracing the function and handle it as a black box.
=#
ModelingToolkit.@register_symbolic PSSE_QUAD_SE(u, SE1, SE2, E1, E2)

#=
On structural prameters:
ideally we'd be able to define the input swicht base one pmech_input and efd_input
as structrual parameters in PSSE_BaseMachine, but that's not possible (don't get properly
forwarde on @extend)
Instead, we define the intermediate varaibles "pmech(t)" and "efd(t)" in base machine
but the parameters get created here and the equation switching also happens here
=#

@mtkmodel PSSE_GENROU begin
    @structural_parameters begin
        pmech_input = true
        efd_input = true
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
        # we cannot @extend a model with structural parameters, so we need to define it in concret mode
        if !pmech_input
            pmech_set, [guess=1, description="mechanical power setpoint [pu]"]
        end
        if !efd_input
            efd_set, [guess=1, description="field voltage setpoint [pu]"]
        end
        # Additional machine parameters for GENROU
        Xpq, [guess=0.6, description="q-axis transient reactance [pu]"]
        Tpq0, [guess=0.9, description="q-axis transient open-circuit time constant [s]"]
    end

    @variables begin
        # State variables (4 additional states for GENROU)
        Epd(t), [guess=0, description="d-axis voltage behind transient reactance [pu]"]
        Epq(t), [guess=1, description="q-axis voltage behind transient reactance [pu]"]
        PSIkd(t), [guess=1, description="d-axis rotor flux linkage [pu]"]
        PSIkq(t), [guess=0, description="q-axis rotor flux linkage [pu]"]

        # Algebraic variables
        PSId(t), [guess=1, description="d-axis flux linkage [pu]"]
        PSIq(t), [guess=0, description="q-axis flux linkage [pu]"]
        PSIppd(t), [guess=1, description="d-axis subtransient flux linkage [pu]"]
        PSIppq(t), [guess=0, description="q-axis subtransient flux linkage [pu]"]
        PSIpp(t), [guess=1, description="Air-gap flux [pu]"]
        XadIfd(t), [guess=1, description="d-axis machine field current [pu]"]
        XaqIlq(t), [guess=0, description="q-axis machine field current [pu]"]
    end

    begin
        # Constants (from OpenIPSL lines 109-116)
        K1d = (Xpd - Xppd)*(Xd - Xpd)/(Xpd - Xl)^2
        K2d = (Xpd - Xl)*(Xppd - Xl)/(Xpd - Xppd)
        K1q = (Xpq - Xppq)*(Xq - Xpq)/(Xpq - Xl)^2
        K2q = (Xpq - Xl)*(Xppq - Xl)/(Xpq - Xppq)
        K3d = (Xppd - Xl)/(Xpd - Xl)
        K4d = (Xpd - Xppd)/(Xpd - Xl)
        K3q = (Xppq - Xl)/(Xpq - Xl)
        K4q = (Xpq - Xppq)/(Xpq - Xl)

        # Alias
        Xpp = Xppd
    end

    @equations begin
        # Input switching (we cannot extend a model with structural parameters)
        pmech ~ pmech_input ? PMECH_in.u : pmech_set
        efd ~ efd_input ? EFD_in.u : efd_set

        # State equations (from OpenIPSL lines 130-133)
        Dt(Epq) ~ 1/Tpd0*(efd - XadIfd)
        Dt(Epd) ~ 1/Tpq0*(-1)*XaqIlq
        Dt(PSIkd) ~ 1/Tppd0*(Epq - PSIkd - (Xpd - Xl)*id)
        Dt(PSIkq) ~ 1/Tppq0*(Epd - PSIkq + (Xpq - Xl)*iq)

        # Override electrical torque calculation (from OpenIPSL line 134)
        Te ~ PSId*iq - PSIq*id

        # Flux linkage equations (from OpenIPSL lines 135-138)
        PSId ~ PSIppd - Xppd*id
        PSIq ~ (-PSIppq) - Xppq*iq
        PSIppd ~ Epq*K3d + PSIkd*K4d
        -PSIppq ~ (-Epd*K3q) - PSIkq*K4q

        # Air-gap flux magnitude (from OpenIPSL line 139)
        PSIpp ~ sqrt(PSIppd*PSIppd + PSIppq*PSIppq)

        # Field current equations with saturation (from OpenIPSL lines 140-151)
        XadIfd ~ K1d*(Epq - PSIkd - (Xpd - Xl)*id) + Epq + id*(Xd - Xpd) +
                 PSSE_QUAD_SE(PSIpp, S10, S12, 1, 1.2)*PSIppd

        XaqIlq ~ K1q*(Epd - PSIkq + (Xpq - Xl)*iq) + Epd - iq*(Xq - Xpq) -
                 PSSE_QUAD_SE(PSIpp, S10, S12, 1, 1.2)*(-1)*PSIppq*(Xq - Xl)/(Xd - Xl)

        # Override voltage equations (from OpenIPSL lines 153-154)
        ud ~ (-PSIq) - R_a*id
        uq ~ PSId - R_a*iq

        # Output connections (from OpenIPSL lines 125-129)
        XADIFD_out.u ~ XadIfd
        ISORCE_out.u ~ XadIfd
    end
end
