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

"""
    PSSE_EXP_SE(u, S_EE_1, S_EE_2, E_1, E_2)

Exponential Saturation Function (PTI PSS/E).
Port of OpenIPSL.NonElectrical.Functions.SE_exp
"""
function PSSE_EXP_SE(u, S_EE_1, S_EE_2, E_1, E_2)
    X = log(S_EE_2/S_EE_1)/log(E_2)
    return S_EE_1*u^X
end

@mtkmodel SimpleLag begin
    @structural_parameters begin
        K # Gain
        T # Time constant
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=0, description="Output signal", output=true]
    end
    @equations begin
        T * Dt(out) ~ K*in - out
    end
end

@mtkmodel SimpleLagLim begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        outMin # Lower limit
        outMax # Upper limi)t
    end
    @parameters begin
        _callback_sat_max
        _callback_sat_min
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=0, description="limited integrator output state"]
        min(t), [description="Lower limit"]
        max(t), [description="Upper limit"]
        forcing(t)
    end
    @equations begin
        min ~ outMin
        max ~ outMax
        forcing ~ K*in - out
        T*Dt(out) ~ (1 - _callback_sat_max - _callback_sat_min) * forcing
    end
end

function attach_SimpleLagLim_callbacks!(cf)
    laglim_components = String[]
    regex = r"^(.*)₊_callback_sat_max$"
    for s in NetworkDynamics.psym(cf)
        m = match(regex, String(s))
        isnothing(m) || push!(laglim_components, m.captures[1])
    end
    if NetworkDynamics.has_callback(cf)
        allcb = NetworkDynamics.get_callbacks(cf)
        for cb in allcb
            if cb isa VectorContinuousComponentCallback && any(s -> !isnothing(match(regex, string(s))), cb.condition.psym)
                error("Component model already has a SimpleLagLim callback attached! Can't attach another one.")
            end
        end
    end
    for ns in laglim_components
        cb = _generate_SimpleLagLim_callbacks(cf, ns)
        NetworkDynamics.add_callback!(cf, cb)
        NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_max"), 0.0)
        NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_min"), 0.0)
    end
    cf
end
function _generate_SimpleLagLim_callbacks(cf::NetworkDynamics.ComponentModel, namespace)
    min = Symbol(namespace, "₊min")
    max = Symbol(namespace, "₊max")
    out = Symbol(namespace, "₊out")
    forcing = Symbol(namespace, "₊forcing")
    satmax = Symbol(namespace, "₊_callback_sat_max")
    satmin = Symbol(namespace, "₊_callback_sat_min")

    condition = ComponentCondition(_SatLim_condition, [min, max, out, forcing], [satmax, satmin])

    upcrossing_affect = ComponentAffect([out], [satmax, satmin]) do u, p, eventidx, ctx
        if eventidx == 1 || eventidx == 2
            println("SimpleLagLim: Enabling max sat at time ", ctx.t)
            p[satmax] = 1.0
            p[satmin] = 0.0
        elseif eventidx == 2
            println("SimpleLagLim: Enabling min sat at time ", ctx.t)
            p[satmax] = 0.0
            p[satmin] = 1.0
        elseif eventidx == 3
            # upcrossing means, forcing went from negative to positive, i.e. we leave lower saturation
            insatmin = !iszero(p[satmin])
            if insatmin
                println("SimpleLagLim: Disabling lower saturation at time ", ctx.t)
                p[satmin] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    downcrossing_affect = ComponentAffect([out],[satmax]) do u, p, eventidx, ctx
        if eventidx == 1 || eventidx == 2
            # in theory should never be hit
            return
        elseif eventidx == 3
            # downcrossing means, forcing went from positive to negative, i.e. we leave upper saturation
            insatmax = !iszero(p[satmax])
            if insatmax
                println("SimpleLagLim: Disabling upper saturation at time ", ctx.t)
                p[satmax] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    discrete_condition = ComponentCondition([out, min, max], []) do u, p, t
        u[out] < u[min] || u[out] > u[max]
    end
    discrete_affect = ComponentAffect([out],[]) do u, p, ctx
        if ctx.model isa VertexModel
            minidx = VIndex(ctx.eidx, min)
            maxidx = VIndex(ctx.eidx, max)
        else
            minidx = EIndex(ctx.eidx, min)
            maxidx = EIndex(ctx.eidx, max)
        end
        min, max = NWState(ctx.integrator)[(minidx, maxidx)]
        if u[out] < min
            @warn "Sanity check cb for LagLim triggered! out=$(u[out]) < min=$min at time $(ctx.t). Forcing out to min. \
                   This might indicate a discrete jump in you model which was not picked up by the callback system!"
            u[out] = min
            p[satmin] = 1.0
            p[satmax] = 0.0
        elseif u[out] > max
            @warn "Sanity check cb for LagLim triggered! out=$(u[out]) > max=$max at time $(ctx.t). Forcing out to max. \
                   This might indicate a discrete jump in you model which was not picked up by the callback system!"
            u[out] = max
            p[satmin] = 0.0
            p[satmax] = 1.0
        else
            error("Sanity check was wrongfully triggered!")
        end
    end

    (
        VectorContinuousComponentCallback(condition, upcrossing_affect, 3; affect_neg! = downcrossing_affect),
        DiscreteComponentCallback(discrete_condition, discrete_affect)
    )
end
function _SatLim_condition(_out, u, p, _)
        # define condition in separate function to avoid capturing and make them batch compatible
        # expect u[min, max, out, forcing]
        # expect p[satmax, satmin]
        insatmax = !iszero(p[1])
        insatmin = !iszero(p[2])

        upcrossing_max =  u[3] - u[2]
        upcrossing_min = -u[3] + u[1]

        # enable upper saturation
        _out[1] = insatmax ? Inf : upcrossing_max
        # enable lower saturation
        _out[2] = insatmin ? Inf : upcrossing_min
        if insatmax || insatmin
            # when in saturation, check for zero crossing of forcing
            # when not in saturation don't set out[3] at all
            _out[3] = u[4]
        end
end

# after Modelica.Blocks.Continuous.Derivative
@mtkmodel Derivative begin
    @structural_parameters begin
        K # Gain
        T # Time constant
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal"]
        internal(t), [guess=0, description="Internal integrator for derivative estimation"]
    end
    @equations begin
        T*Dt(internal) ~ in - internal
        out ~ K/T*(in - internal)
    end
end

