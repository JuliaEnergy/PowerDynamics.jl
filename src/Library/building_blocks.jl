"""
    QUAD_SE(u, SE1, SE2, E1, E2)

Quadratic Saturation Function through two points (E1,SE1) and (E2,SE2).
"""
function QUAD_SE(u, SE1, SE2, E1, E2)
    if any(x -> x isa GradientTracer, (u, SE1, SE2, E1, E2))
        return u+SE1+SE2+E1+E2
    end
    # the check is annoying the model might get tested with arbitrary data...
    # if !(0 < SE1 < SE2) || !(0 < E1 < E2)
    #     throw(ArgumentError("QUAD_SE: Saturation values and voltage points must be positive and increasing! Got SE1=$SE1, SE2=$SE2, E1=$E1, E2=$E2"))
    # end
    SE1 == SE2 && return SE1

    a = sqrt(SE1 * E1 / (SE2 * E2))
    A = E2 - (E1 - E2) / (a - 1)
    if u <= A
        return 0.0
    else
        B = SE2 * E2 * (a - 1)^2 / (E1 - E2)^2
        return B * (u - A)^2 / u
    end
end
ModelingToolkit.@register_symbolic QUAD_SE(u, SE1, SE2, E1, E2)

"""
    EXP_SE(u, SE1, SE2, E1, E2)

Exponential Saturation Function through two points (E1,SE1) and (E2,SE2).
"""
function EXP_SE(u, SE1, SE2, E1, E2)
    if any(x -> x isa GradientTracer, (u, SE1, SE2, E1, E2))
        return u+SE1+SE2+E1+E2
    end
    # commented out check, to not register as symbolic
    # if !(0 < SE1 < SE2) || !(0 < E1 < E2)
    #     throw(ArgumentError("EXP_SE: Saturation values and voltage points must be positive and increasing! Got SE1=$SE1, SE2=$SE2, E1=$E1, E2=$E2"))
    # end
    SE1 == SE2 && return SE1

    X = log(SE2/SE1) / log(E2/E1)
    k = SE1 / E1^X

    return k * u^X  # equivalently: SE1 * (u/E1)^X
end
ModelingToolkit.@register_symbolic EXP_SE(u, SE1, SE2, E1, E2)

"""
SimpleLag block

```asciiart
    ╭─────────╮
 in │    K    │ out
╶───┤╶───────╴├────╴
    │ 1 + s T │
    ╰─────────╯
```
Additional structural parameters:
- `guess=0`/`default`: initial guess/default for the internal state (equals output in steady state)
- `allowzeroT`: if true, the lag is be bypassed when T=0 (this does not reduce the model order)
"""
@mtkmodel SimpleLag begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        guess=0
        default=nothing
        allowzeroT=false
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        internal(t)=default, [guess=guess]
        out(t), [description="Output signal", output=true]
    end
    @equations begin
        if allowzeroT
            Dt(internal) ~ ifelse(T==0, 0, (K*in - internal)/T)
            out ~ ifelse(T==0, in, internal)
        else
            T * Dt(internal) ~ K*in - internal
            out ~ internal
        end
    end
end

"""
SimpleLead block

```asciiart
    ╭─────────╮
 in │ 1 + s T │ out
╶───┤╶───────╴├────╴
    │    K    │
    ╰─────────╯
```

This block directly uses `Dt(in)`, therefore it does not add additional states
but may not be used in all scenarios!
"""
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


"""
    SaturationConfig(method; discrete_cb=true, continuous_cb=true, regularization=0)

Configuration for saturation limiting and anti-windup in limited integrator blocks
([`SimpleLagLim`](@ref), [`LimIntegrator`](@ref)).

This controls how output limits are enforced and how integrator windup is prevented when
limits are active. Different methods have different numerical properties and may affect
initialization and simulation accuracy.

# Supported methods
- `:callback`: Use callbacks to detect and enforce limits, stopping integration when saturated (default)
- `:complementary`: Use complementarity constraints (Fischer-Burmeister) with Lagrange multipliers
- `:rhs_hard`: Clamp RHS with hard `ifelse` switching (discontinuous anti-windup)
- `:rhs_soft`: Clamp RHS with smooth tanh-based saturation (smooth anti-windup)

# Arguments
- `method`: Saturation handling and anti-windup method
- `discrete_cb=true`: Enable discrete callbacks (only for `:callback` method)
- `continuous_cb=true`: Enable continuous callbacks (only for `:callback` method)
- `regularization=0`: Regularization parameter (for `:complementary` and `:rhs_soft` methods)

See also [`SaturationConfiguration`](@ref), [`set_saturation_config!`](@ref).
"""
struct SaturationConfig
    method::Symbol
    discrete_cb::Bool
    continuous_cb::Bool
    regularization::Float64
end
function SaturationConfig(method; discrete_cb=true, continuous_cb=true, regularization=0, verbose=true)
    if method == :callback
        @assert regularization == 0 "Regularization not supported for method :callback"
        return SaturationConfig(method, discrete_cb, continuous_cb, regularization)
    elseif method == :complementary
        return SaturationConfig(method, false, false, regularization)
    elseif method == :rhs_hard
        return SaturationConfig(method, false, false, NaN)
    elseif method == :rhs_soft
        return SaturationConfig(method, false, false, regularization)
    else
        error("Unknown saturation config method :$method. Supported types are :callback, :complementary, :rhs_hard and :rhs_soft")
    end
end

"""
    const SaturationConfiguration = ScopedValue(SaturationConfig(:callback))

Global configuration for saturation handling in limited integrator blocks.
See [`SaturationConfig`](@ref) for details on available options.

Use [`set_saturation_config!`](@ref) to change globally or use
```julia
with(SaturationConfiguration => SaturationConfig(:complementary, regularization=1e-6)) do
    # your code here
end
```
to set temporarily via ScopedValue mechanism.
"""
const SaturationConfiguration = ScopedValue(SaturationConfig(:callback))

"""
    set_saturation_config!(style::Symbol; kwargs...)
    set_saturation_config!(config::SaturationConfig)

Sets [`SaturationConfiguration`](@ref) globally.
See [`SaturationConfig`](@ref) for available options and styles.
"""
function set_saturation_config!(config::SaturationConfig)
    SaturationConfiguration[] = config
end
function set_saturation_config!(style::Symbol; kwargs...)
    config = SaturationConfig(style; kwargs...)
    set_saturation_config!(config)
end

@component function LimitedIntegratorBase(; name, type, K, T, outMin, outMax, guess=0)
    if type != :lag && type != :int
        error("Unknown type $type for SimpleLagLim. Supported types are :lag and :int")
    end
    config = SaturationConfiguration[]
    if config.method ∉ (:callback, :complementary, :rhs_hard, :rhs_soft)
        error("Unknown saturation config method :$(config.method). Supported types are :callback, :complementary, :rhs_hard and :rhs_soft")
    end

    @parameters begin
        # only used for callbacks
        _callback_sat_max, [guess=0, description="internal callback parameter indicating upper saturation state"]
        _callback_sat_min, [guess=0, description="internal callback parameter indicating lower saturation state"]
        # onlye used for rhs_soft and complementary
        ϵ=config.regularization, [description="Regularization parameter for saturation handling"]
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        x(t), [guess=guess, description="limited integrator output state"]
        out(t), [description="limited integrator output", output=true]
        min(t), [description="Lower limit"]
        max(t), [description="Upper limit"]
        forcing(t), [description="rhs of internal integrator"]
        # only used in case of complementary constraints
        λ⁻(t), [guess=0, description="lagrange multiplier for lower limit", bounds=(0, Inf)]
        λ⁺(t), [guess=0, description="lagrange multiplier for upper limit", bounds=(0, Inf)]
    end

    forcing_eqs = if type == :lag
        forcing ~ K*in - x
    elseif type == :int
        forcing ~ K*in
    end

    eqs = [
        min ~ outMin
        max ~ outMax
        forcing_eqs
    ]

    # integrator eqations T*D(x) ~ forcing depend on config method
    if config.method == :callback
        neweqs = [
            T*Dt(x) ~ (1 - _callback_sat_max - _callback_sat_min) * forcing
            out ~ x
        ]
    elseif config.method == :complementary
        # Fischer–Burmeister constraints
        eqλ⁻ = 0 ~ sqrt(λ⁻^2 + (x - min)^2 + ϵ) - λ⁻ - (x - min)
        eqλ⁺ = 0 ~ sqrt(λ⁺^2 + (max - x)^2 + ϵ) - λ⁺ - (max - x)
        # Dynamics with multipliers depending on which limits are active
        neweqs = if !_isneginf(outMin)  && !_isposinf(outMax)
            [
                T * Dt(x) ~ forcing - λ⁺ + λ⁻
                eqλ⁻
                eqλ⁺
            ]
        elseif _isneginf(outMin) && !_isposinf(outMax)
            [
                T * Dt(x) ~ forcing - λ⁺
                eqλ⁺
            ]
        elseif !_isneginf(outMin)  && _isposinf(outMax)
            [
                T * Dt(x) ~ forcing + λ⁻
                eqλ⁻
            ]
        elseif _isneginf(outMin) && _isposinf(outMax)
            [
                T * Dt(x) ~ forcing
            ]
        end
        push!(neweqs, out ~ x)
    elseif config.method == :rhs_hard
        neweqs = [
            T*Dt(x) ~  _hard_clamped_rhs(forcing, x, min, max)
            out ~ clamp(x, min, max)
        ]
    elseif config.method == :rhs_soft
        neweqs = [
            T*Dt(x) ~  _soft_clamped_rhs(forcing, x, min, max, ϵ)
            out ~ clamp(x, min, max)
        ]
    end
    append!(eqs, neweqs)

    sys = System(eqs, t; name)

    sys = if config.method == :callback
        setmetadata(sys, ComponentPostprocessing, attach_limint_postprocessing_callback!)
    elseif config.method == :complementary
        setmetadata(sys, ComponentPostprocessing, attach_limint_postprocessing_complementary!)
    elseif config.method ∈ (:rhs_hard, :rhs_soft)
        setmetadata(sys, ComponentPostprocessing, attach_limint_postprocessing_rhs!)
    else
        sys
    end
    return sys
end
_isposinf(x::Num) = false
_isneginf(x::Num) = false
_isposinf(x::Float64) = x == Inf
_isneginf(x::Float64) = x == -Inf
_isposinf(x::Int) = x == typemax(Int)
_isneginf(x::Int) = x == -typemin(Int)
function _soft_clamped_rhs(u, x, xmin, xmax, ε)
    s_hi = 0.5 * (1 - tanh((x - xmax)/ε))
    s_lo = 0.5 * (1 + tanh((x - xmin)/ε))
    # Only scale positive forcing by s_hi (upper limit)
    # Only scale negative forcing by s_lo (lower limit)
    return ifelse(u > 0, u * s_hi, u * s_lo)
end
function _hard_clamped_rhs(u, x, xmin, xmax)
    # make compatible with sparacity tracing
    if any(x -> x isa GradientTracer, (u, x, xmin, xmax))
        return u+x+xmin+xmax
    end
    ifelse(((x ≥ xmax) && (u > 0)) || ((x ≤ xmin) && (u < 0)), 0.0, u)
end
Symbolics.@register_symbolic _hard_clamped_rhs(u, x, xmin, xmax)

function attach_limint_postprocessing_callback!(cf, ns)
    cb = _generate_limint_callbacks(ns)
    NetworkDynamics.add_callback!(cf, cb)
    NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_max"), 0.0)
    NetworkDynamics.set_default!(cf, Symbol(ns, "₊_callback_sat_min"), 0.0)

    #=
    # TODO: Experiment: add constraints and out ~ clamp(x, min, max) here to enforce consistency at initialization
    # maybe look into Juniper and solve als MINLP problem?
    x = Symbol(ns, :₊x)
    out = Symbol(ns, :₊out)
    min = Symbol(ns, :₊min)
    max = Symbol(ns, :₊max)
    satmin = Symbol(ns, :₊_callback_sat_min)
    satmax = Symbol(ns, :₊_callback_sat_max)

    # Force x to be within bounds: combined with out=clamp(x,min,max), this ensures x∈[min,max]
    # Saturation flags must be 0 or 1 to properly stop dynamics when saturated
    forcing = Symbol(ns, :₊forcing)
    ic = NetworkDynamics.InitConstraint([x, out, min, max, forcing, satmin, satmax], 2) do res, u
        res[1] = u[x] - clamp(u[x], u[min], u[max])
        # res[2] = u[satmax] - (u[x] >= u[max] > 0 ? 1.0 : 0.0)
        # res[3] = u[satmin] - (u[x] <= u[min] < 0 ? 1.0 : 0.0)
        res[2] = sqrt(u[satmin]^2 + u[satmax]^2 + 1e-6) - u[satmin] - u[satmax]
        nothing
    end
    NetworkDynamics.add_initconstraint!(cf, ic)
    =#
end
function attach_limint_postprocessing_complementary!(cf, ns)
    # # add explicit limits
    # soft_min(a, b, δ) = -δ * log(exp(-a/δ) + exp(-b/δ))

    # λ⁺ = Symbol(ns, :₊λ⁺)
    # λ⁻ = Symbol(ns, :₊λ⁻)
    # min = Symbol(ns, :₊min)
    # max = Symbol(ns, :₊max)
    # x = Symbol(ns, :₊x)
    # ϵ = Symbol(ns, :₊ϵ)

    # if λ⁺ ∈ NetworkDynamics.sym(cf) # has upper bound
    #     upperc = NetworkDynamics.InitConstraint([x, max, ϵ], 1) do res, u
    #         # Logarithmic barrier: penalizes violations strongly
    #         # Goes to +∞ as x approaches max from above
    #         # res[1] = Base.min(0.0, u[max] - u[x])
    #         res[1] = -log(u[max] - u[x] + sqrt(u[ϵ]))
    #     end
    #     NetworkDynamics.add_initconstraint!(cf, upperc)
    # end
    # if λ⁻ ∈ NetworkDynamics.sym(cf) # has lower bound
    #     lowerc = NetworkDynamics.InitConstraint([x, min, ϵ], 1) do res, u
    #         # Logarithmic barrier: penalizes violations strongly
    #         # Goes to +∞ as x approaches min from below
    #         # res[1] = Base.min(0.0, u[x] - u[min])
    #         res[1] = -log(u[x] - u[min] + sqrt(u[ϵ]))
    #     end
    #     NetworkDynamics.add_initconstraint!(cf, lowerc)
    # end
end
function attach_limint_postprocessing_rhs!(cf, ns)
    x = Symbol(ns, :₊x)
    out = Symbol(ns, :₊out)
    ic = NetworkDynamics.InitConstraint([x, out], 1) do res, u
        res[1] = u[out] - u[x]
    end
    NetworkDynamics.add_initconstraint!(cf, ic)
end

"""
SimpleLagLim block

```asciiart
              __ outMax
             /
      ╭─────────╮
   in │    K    │ out
  ╶───┤╶───────╴├────╴
      │ 1 + s T │
      ╰─────────╯
outMin __/
```

Additional structural parameters:
- `guess=0`: initial guess for the internal state (equals output in steady state)
"""
SimpleLagLim(; kwargs...) = LimitedIntegratorBase(; type=:lag, kwargs...)

"""
LimIntegrator block

```asciiart
              __ outMax
             /
        ╭─────╮
     in │  K  │ out
    ╶───┤╶───╴├────╴
        │ s T │
        ╰─────╯
outMin __/
```

Additional structural parameters:
- `guess=0`: initial guess for the internal state (equals output in steady state)
"""
LimIntegrator(; kwargs...) = LimitedIntegratorBase(; type=:int, T=1, kwargs...)

function _generate_limint_callbacks(namespace)
    min = Symbol(namespace, "₊min")
    max = Symbol(namespace, "₊max")
    x = Symbol(namespace, "₊x")
    forcing = Symbol(namespace, "₊forcing")
    satmax = Symbol(namespace, "₊_callback_sat_max")
    satmin = Symbol(namespace, "₊_callback_sat_min")

    use_continuous_callback = SaturationConfiguration[].continuous_cb
    use_discrete_callback = SaturationConfiguration[].discrete_cb

    if use_continuous_callback
        condition = ComponentCondition(_SatLim_condition, [min, max, x, forcing], [satmax, satmin])

        upcrossing_affect = ComponentAffect([], [satmax, satmin]) do u, p, eventidx, ctx
            comp = get_compidx(ctx)
            verbose = CallbackVerbose[]
            if eventidx == 1
                verbose && println("$comp: $namespace: /⎺ reached upper saturation at $(round(ctx.t, digits=4))s")
                p[satmax] = 1.0
                p[satmin] = 0.0
            elseif eventidx == 2
                verbose && println("$comp: $namespace: \\_ reached lower saturation at $(round(ctx.t, digits=4))s")
                p[satmax] = 0.0
                p[satmin] = 1.0
            elseif eventidx == 3
                # upcrossing means, forcing went from negative to positive, i.e. we leave lower saturation
                insatmin = !iszero(p[satmin])
                if insatmin
                    verbose && println("$comp: $namespace: _/ left lower saturation at $(round(ctx.t, digits=4))s")
                    p[satmin] = 0.0
                end
            else
                error("Unknown event index $eventidx")
            end
        end

        downcrossing_affect = ComponentAffect([],[satmax]) do u, p, eventidx, ctx
            comp = get_compidx(ctx)
            verbose = CallbackVerbose[]
            if eventidx == 1 || eventidx == 2
                # in theory should never be hit
                return
            elseif eventidx == 3
                # downcrossing means, forcing went from positive to negative, i.e. we leave upper saturation
                insatmax = !iszero(p[satmax])
                if insatmax
                    verbose && println("$comp: $namespace: ⎺\\ left upper saturation at $(round(ctx.t, digits=4))s")
                    p[satmax] = 0.0
                end
            else
                error("Unknown event index $eventidx")
            end
        end
        continous_callback = VectorContinuousComponentCallback(condition, upcrossing_affect, 3; affect_neg! = downcrossing_affect)
    end

    if use_discrete_callback
        # TODO: merge both discrete conditions and move condition below function for performance
        function _discrete_cond(u,p,t)
            # account for nummerical innaccuracies at the boudaries
            u[1] < u[2] - 1e-10 || u[1] > u[3] + 1e-10
        end
        discrete_condition = ComponentCondition(_discrete_cond, [x, min, max], [])
        discrete_affect = ComponentAffect([x],[satmin, satmax]) do u, p, ctx
            comp = get_compidx(ctx)
            verbose = CallbackVerbose[]
            if ctx.model isa VertexModel
                minidx = VIndex(ctx.vidx, min)
                maxidx = VIndex(ctx.vidx, max)
            else
                minidx = EIndex(ctx.eidx, min)
                maxidx = EIndex(ctx.eidx, max)
            end
            _min, _max = NWState(ctx.integrator)[(minidx, maxidx)]
            if u[x] < _min
                if verbose || use_continuous_callback
                    print("$comp: $namespace: \\_ reached lower saturation at $(round(ctx.t, digits=4))s")
                    use_continuous_callback ? printstyled(" (triggered by discrete cb)\n", color=:yellow) : println()
                end
                u[x] = _min
                p[satmin] = 1.0
                p[satmax] = 0.0
            elseif u[x] > _max
                if verbose || use_continuous_callback
                    print("$comp: $namespace: /⎺ reached upper saturation at $(round(ctx.t, digits=4))s")
                    use_continuous_callback ? printstyled(" (triggered by discrete cb)\n", color=:yellow) : println()
                end
                u[x] = _max
                p[satmin] = 0.0
                p[satmax] = 1.0
            else
                error("Sanity check was wrongfully triggered!")
            end
        end
        function _discrete_unsat_cond(u,p,t)
            insatmin = !iszero(p[1])
            insatmax = !iszero(p[2])
            insatmin && u[1] > 0 || insatmax && u[1] < 0
        end
        discrete_unsat_condition = ComponentCondition(_discrete_unsat_cond, [forcing],[satmin, satmax])
        discrete_unsat_affect = ComponentAffect([],[satmin, satmax]) do u, p, ctx
            comp = get_compidx(ctx)
            verbose = CallbackVerbose[]
            insatmin = !iszero(p[satmin])
            insatmax = !iszero(p[satmax])
            if insatmin
                if verbose
                    print("$comp: $namespace: _/ left lower saturation at $(round(ctx.t, digits=4))s")
                    use_continuous_callback ? printstyled(" (triggered by discrete cb)\n", color=:yellow) : println()
                end
                p[satmin] = 0.0
            elseif insatmax
                if verbose
                    print("$comp: $namespace: ⎺\\ left upper saturation at $(round(ctx.t, digits=4))s")
                    use_continuous_callback ? printstyled(" (triggered by discrete cb)\n", color=:yellow) : println()
                end
                p[satmax] = 0.0
            else
                error("Sanity check was wrongfully triggered!")
            end
        end

        discrete_callback = (
            DiscreteComponentCallback(discrete_condition, discrete_affect),
            DiscreteComponentCallback(discrete_unsat_condition, discrete_unsat_affect)
        )
    end

    if use_continuous_callback && use_discrete_callback
        return (continous_callback, discrete_callback...)
    elseif use_continuous_callback
        return continous_callback
    elseif use_discrete_callback
        return discrete_callback
    else
        error("No callback type selected for LimiterCallbackConfig!")
    end
end
function _SatLim_condition(_out, u, p, _)
        # define condition in separate function to avoid capturing and make them batch compatible
        # expect u[min, max, x, forcing]
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
            _out[3] = u[4]
        else
            # when not in saturation, set out[3] at Inf
            # This might be problematic if
            # - lower lim is hit (i.e. forcing is negative)
            # - next round, forcing is still negativ so we have a discrete jump from Inf to small negativ, which is a zero crossing
            # - but it seems like this non-continuous crossing is not actually registered as a crossing? Maybe because t=t in both cases?
            _out[3] = Inf
        end
end
get_compidx(nt::NamedTuple) = haskey(nt, :vidx) ? VIndex(nt.vidx) : EIndex(nt.eidx)

"""
Simple gain block

```asciiart
 in ╭───╮ out
╶───┤ K ├────╴
    ╰───╯
```
"""
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

"""
Derivative approximation block. Modeld after Modelica.Blocks.Continuous.Derivative

```asciiart
    ╭─────────╮
 in │   s K   │ out
╶───┤╶───────╴├────╴
    │ 1 + s T │
    ╰─────────╯
```

Additional structural parameters:
- `guess=0`: initial guess for the internal state (equals input in steady state)
"""
@mtkmodel Derivative begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        guess=0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal"]
        internal(t), [guess=guess, description="Internal integrator for derivative estimation"]
    end
    @equations begin
        T*Dt(internal) ~ in - internal
        out ~ K/T*(in - internal)
    end
end

"""
LeadLag block

```asciiart
    ╭──────────╮
 in │  1 + sT1 │ out
╶───┤K╶───────╴├────╴
    │  1 + sT2 │
    ╰──────────╯
```

Additional structural parameters:
- `guess=0`: initial guess for the internal state (equals input in steady state)
- `allowzeroT`: if true, the lead-lag is be bypassed when T1=0 and T2=0 (this does not reduce the model order)
"""
@mtkmodel LeadLag begin
    @structural_parameters begin
        K # Gain
        T1 # Lead time constant
        T2 # Lag time constant
        guess=0
        allowzeroT=false
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [description="Output signal", output=true]
        internal(t), [guess=guess, description="Internal state"]
        internal_dt(t), [description="derivative of internal state"]
    end
    @equations begin
        if allowzeroT
            internal_dt ~ ifelse((T1==0) & (T2==0), 0, (in - internal)/T2)
            out ~ ifelse((T1==0) & (T2==0), K*in, K*(internal + T1*internal_dt))
        else
            internal_dt ~ (in - internal)/T2
            out ~ K*(internal + T1*internal_dt)
        end
        Dt(internal) ~ internal_dt
    end
end

"""
DeadZone block, modeled after Modelica.Blocks.Nonlinear.DeadZone

```asciiart
        │    ╱
   uMin │   ╱
─────┼╼━┿━╾┼─────
    ╱   │ uMax
   ╱    │
```

A dead zone nonlinearity that outputs zero when the input is within the specified band [uMin, uMax].
- If `in < uMin`: `out = in - uMin` (negative linear)
- If `uMin ≤ in ≤ uMax`: `out = 0` (dead zone)
- If `in > uMax`: `out = in - uMax` (positive linear)

Structural parameters:
- `uMax`: Upper dead zone limit
- `uMin=-uMax`: Lower dead zone limit (defaults to -uMax for symmetric dead zone)
"""
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

"""
    ss_to_mtkmodel(A, B, C, D; name=nothing, guesses=zeros(size(A,1)))

Convert a state-space representation to a ModelingToolkit model.

Matrices can be either of real numbers or symbolic parameters/terms.

Returns A `System` object with variables `in` (input), `out` (output), and `x₁, x₂, ...` (states).
"""
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

"""
    siso_tf_to_ss(num, den)

Convert a SISO transfer function to state-space representation.

Takes polynomial coefficients for numerator and denominator and returns state-space matrices (A, B, C, D).
The transfer function is represented as:
```
       num[1]sⁿ + num[2]sⁿ⁻¹ + ... + num[end]
G(s) = ─────────────────────────────────────
       den[1]sᵐ + den[2]sᵐ⁻¹ + ... + den[end]
```

# Arguments
- `num`: Vector of numerator coefficients (highest degree first)
- `den`: Vector of denominator coefficients (highest degree first)

# Returns
A tuple `(A, B, C, D)` of state-space matrices in controller canonical form.

# Notes
- Leading zeros in num/den are automatically truncated
- Transfer function must be proper (numerator degree ≤ denominator degree)
- Adapted from SymbolicControlSystems.jl (Copyright (c) 2020 Fredrik Bagge Carlson, MIT License)
"""
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
        C[:] .-= bN .* reverse(den)[1:end-1] # Can index into polynomials at greater indices than their length
    end
    D = fill(bN, 1, 1)

    return A, B, C, D
end
