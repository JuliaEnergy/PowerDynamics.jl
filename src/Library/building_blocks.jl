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
        guess=0
        default=nothing
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t)=default, [guess=guess, description="Output signal", output=true]
    end
    @equations begin
        T * Dt(out) ~ K*in - out
    end
end
@mtkmodel SimpleLead begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        guess=0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=guess, description="Output signal", output=true]
    end
    @equations begin
        T*Dt(in) ~ K*out - in
    end
end

@mtkmodel LimitedIntegratorBase begin
    @structural_parameters begin
        type # :lag or :int
        K # Gain
        T # Time constant
        outMin # Lower limit
        outMax # Upper limit
        guess=0
    end
    @parameters begin
        _callback_sat_max
        _callback_sat_min
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=guess, description="limited integrator output state", output=true]
        min(t), [description="Lower limit"]
        max(t), [description="Upper limit"]
        forcing(t)
    end
    @equations begin
        min ~ outMin
        max ~ outMax
        if type == :lag
            forcing ~ K*in - out
        elseif type == :int
            forcing ~ K*in
        else
            error("Unknown type $type for SimpleLagLim. Supported types are :lag and :int")
        end
        T*Dt(out) ~ (1 - _callback_sat_max - _callback_sat_min) * forcing
    end
end
SimpleLagLim(; kwargs...) = LimitedIntegratorBase(; type=:lag, kwargs...)
LimIntegrator(; kwargs...) = LimitedIntegratorBase(; type=:int, T=1, kwargs...)

function attach_limint_callbacks!(cf)
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
        cb = _generate_limint_callbacks(cf, ns)
        NetworkDynamics.add_callback!(cf, cb)
        NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_max"), 0.0)
        NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_min"), 0.0)
    end
    cf
end
function _generate_limint_callbacks(cf::NetworkDynamics.ComponentModel, namespace)
    min = Symbol(namespace, "₊min")
    max = Symbol(namespace, "₊max")
    out = Symbol(namespace, "₊out")
    forcing = Symbol(namespace, "₊forcing")
    satmax = Symbol(namespace, "₊_callback_sat_max")
    satmin = Symbol(namespace, "₊_callback_sat_min")

    condition = ComponentCondition(_SatLim_condition, [min, max, out, forcing], [satmax, satmin])

    upcrossing_affect = ComponentAffect([], [satmax, satmin]) do u, p, eventidx, ctx
        if eventidx == 1
            println("$namespace: /⎺ reached upper saturation at $(round(ctx.t, digits=4))s")
            p[satmax] = 1.0
            p[satmin] = 0.0
        elseif eventidx == 2
            println("$namespace: \\_ reached lower saturation at $(round(ctx.t, digits=4))s")
            p[satmax] = 0.0
            p[satmin] = 1.0
        elseif eventidx == 3
            # upcrossing means, forcing went from negative to positive, i.e. we leave lower saturation
            insatmin = !iszero(p[satmin])
            if insatmin
                println("$namespace: _/ left lower saturation at $(round(ctx.t, digits=4))s")
                p[satmin] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    downcrossing_affect = ComponentAffect([],[satmax]) do u, p, eventidx, ctx
        if eventidx == 1 || eventidx == 2
            # in theory should never be hit
            return
        elseif eventidx == 3
            # downcrossing means, forcing went from positive to negative, i.e. we leave upper saturation
            insatmax = !iszero(p[satmax])
            if insatmax
                println("$namespace: ⎺\\ left upper saturation at $(round(ctx.t, digits=4))s")
                p[satmax] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    # TODO: merge both discrete conditions and move condition below function for performance
    discrete_condition = ComponentCondition([out, min, max], []) do u, p, t
        # account for nummerical innaccuracies at the boudaries
        u[out] < u[min] - 1e-10 || u[out] > u[max] + 1e-10
    end
    discrete_affect = ComponentAffect([out],[satmin, satmax]) do u, p, ctx
        if ctx.model isa VertexModel
            minidx = VIndex(ctx.vidx, min)
            maxidx = VIndex(ctx.vidx, max)
        else
            minidx = EIndex(ctx.eidx, min)
            maxidx = EIndex(ctx.eidx, max)
        end
        _min, _max = NWState(ctx.integrator)[(minidx, maxidx)]
        if u[out] < _min
            @warn "Sanity check cb for LagLim triggered! out=$(u[out]) < min=$_min at time $(ctx.t). Forcing out to min. \
                   This might indicate a discrete jump in you model which was not picked up by the callback system!"
            u[out] = _min
            p[satmin] = 1.0
            p[satmax] = 0.0
        elseif u[out] > _max
            @warn "Sanity check cb for LagLim triggered! out=$(u[out]) > max=$_max at time $(ctx.t). Forcing out to max. \
                   This might indicate a discrete jump in you model which was not picked up by the callback system!"
            u[out] = _max
            p[satmin] = 0.0
            p[satmax] = 1.0
        else
            error("Sanity check was wrongfully triggered!")
        end
    end
    discrete_unsat_condition = ComponentCondition([forcing],[satmin, satmax]) do u, p, t
        insatmin = !iszero(p[satmin])
        insatmax = !iszero(p[satmax])
        insatmin && u[forcing] > 0 || insatmax && u[forcing] < 0
    end
    discrete_unsat_affect = ComponentAffect([],[satmin, satmax]) do u, p, ctx
        insatmin = !iszero(p[satmin])
        insatmax = !iszero(p[satmax])
        if insatmin
            println("$namespace: _/ left lower saturation at $(round(ctx.t, digits=4))s (triggered by discrete cb)")
            p[satmin] = 0.0
        elseif insatmax
            println("$namespace: ⎺\\ left upper saturation at $(round(ctx.t, digits=4))s (triggered by discrete cb)")
            p[satmax] = 0.0
        else
            error("Sanity check was wrongfully triggered!")
        end
    end

    (
        VectorContinuousComponentCallback(condition, upcrossing_affect, 3; affect_neg! = downcrossing_affect),
        DiscreteComponentCallback(discrete_condition, discrete_affect),
        DiscreteComponentCallback(discrete_unsat_condition, discrete_unsat_affect)
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

@mtkmodel SimpleGain begin
    @structural_parameters begin
        K # Gain
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal", output=true]
    end
    @equations begin
        out ~ K*in
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

# after OpenIPSL.NonElectrical.Continuous.LeadLag
@mtkmodel LeadLag begin
    @structural_parameters begin
        K # Gain
        T1 # Lead time constant
        T2 # Lag time constant
        guess=0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal", output=true]
        internal(t), [guess=guess, description="Internal state"]
        internal_dt(t), [description="derivative of internal state"]
    end
    @equations begin
        internal_dt ~ (in - internal)/T2
        Dt(internal) ~ internal_dt
        out ~ K*(internal + T1*internal_dt)
    end
end

# after modelica Modelica.Blocks.Nonlinear.DeadZone
@mtkmodel DeadZone begin
    @structural_parameters begin
        uMax # Lower dead zone limit
        uMin=-uMax # Upper dead zone limit
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal", output=true]
    end
    @equations begin
        out ~ ifelse(in < uMin,
            in - uMin,
            ifelse(in > uMax,
                in - uMax,
                0.0
            )
        )
    end
end

ModelingToolkit.@component ss_to_mtkmodel(; A, B, C, D, kwargs...) = ss_to_mtkmodel(A, B, C, D; kwargs...)
function ss_to_mtkmodel(A, B, C, D; name=nothing, guesses=zeros(size(A,1)))
    t = ModelingToolkit.t_nounits
    Dt = ModelingToolkit.D_nounits

    n = size(A, 1)
    @assert size(D) == (1, 1) "Only SISO systems supported"

    # Symbolic system
    @variables in(t) out(t)
    _xs_names = [ Symbol("x", NetworkDynamics.subscript(i)) for i in 1:n]
    # dont use Symbolics.variables as it does not create all necessary metadata?
    # also needs to set guess in @variables not with MTK.setguess
    x = map(zip(_xs_names, guesses)) do (_name, _guess)
        only(@variables $(_name)(t) [guess=_guess])
    end

    ∂x = Dt.(x)
    eqs = vcat(
        ∂x .~ A*x .+ B*[in],
        [out] .~ (length(C)>0 ? C*x : 0) .+ D*[in]
    )
    eqs = Symbolics.simplify.(eqs)
    allp = mapreduce(Symbolics.get_variables, ∪, Iterators.flatten((A,B,C,D)))

    return System(eqs, t, vcat(x, [in, out]), allp; name=name)
end

#=
Taken and adapted from SymbolicControlSystems.jl

Copyright (c) 2020 Fredrik Bagge Carlson, MIT License
=#
function siso_tf_to_ss(num0, den0)
    T = Base.promote_type(eltype(num0), eltype(den0))

    # truncate leading zeros
    num0 = num0[findfirst(!iszero, num0):end]
    den0 = den0[findfirst(!iszero, den0):end]

    # check if it is proper
    denorder = length(den0) - 1
    numorder = length(num0) - 1
    if numorder > denorder
        error("Numerator degree > denominator degree not allowed (non-proper).")
    end


    # Normalize the numerator and denominator to allow realization of transfer functions
    # that are proper, but not strictly proper
    num = num0 ./ den0[1]
    den = den0 ./ den0[1]

    N = length(den) - 1 # The order of the rational function f

    # Get numerator coefficient of the same order as the denominator
    bN = length(num) == N+1 ? num[1] : zero(eltype(num))

    @views if N == 0 #|| num == zero(Polynomial{T})
        A = zeros(T, 0, 0)
        B = zeros(T, 0, 1)
        C = zeros(T, 1, 0)
    else
        A = LinearAlgebra.diagm(1 => ones(T, N-1))
        A[end, :] .= .-reverse(den)[1:end-1]

        B = zeros(T, N, 1)
        B[end] = one(T)

        C = zeros(T, 1, N)
        C[1:min(N, length(num))] = reverse(num)[1:min(N, length(num))]
        C[:] .-= bN .* reverse(den)[1:end-1] # Can index into polynomials at greater inddices than their length
    end
    D = fill(bN, 1, 1)

    return A, B, C, D
end
