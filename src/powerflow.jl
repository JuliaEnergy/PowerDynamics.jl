"""
    pfSlack(; V=missing, δ=missing, u_r=missing, u_i=missing, name=:slackbus)

Create a slack bus for power flow analysis.

A slack bus maintains constant voltage magnitude and phase angle (or real and imaginary voltage components).
Either provide voltage magnitude `V` and phase angle `δ`, or provide real and imaginary voltage components `u_r` and `u_i`.
"""
function pfSlack(; V=missing, δ=missing, u_r=missing, u_i=missing, name=:slackbus, kwargs...)
    if !ismissing(V) && ismissing(u_r) && ismissing(u_i)
        δ = ismissing(δ) ? 0 : δ
        @named slack = Library.VδConstraint(; V, δ)
        u_r = V * cos(δ)
        u_i = V * sin(δ)
    elseif ismissing(V) && ismissing(δ) && !ismissing(u_r) && !ismissing(u_i)
        @named slack = Library.UrUiConstraint(; u_r, u_i)
    else
        throw(ArgumentError("Either Provide V or δ, or u_r and u_i. But not both!"))
    end

    mtkbus = MTKBus(slack; name)
    b = compile_bus(mtkbus; kwargs...)
    set_voltage!(b, u_r + im * u_i)
    set_current!(b, -1.0) # probably injection
    b
end

"""
    pfPV(; P, V, name=:pvbus)

Create a PV bus for power flow analysis.

A PV bus maintains constant active power injection and voltage magnitude.
The reactive power and voltage phase angle are determined by the power flow solution.
"""
function pfPV(; P, V, name=:pvbus, kwargs...)
    @named pv = Library.PVConstraint(; P, V)
    mtkbus = MTKBus(pv; name)
    b = compile_bus(mtkbus; kwargs...)
    set_voltage!(b; mag=V, arg=0)
    i = P/V
    set_current!(b, -i)
    set_default!(b, :pv₊terminal₊i_r, i) # necessary for current source
    set_default!(b, :pv₊terminal₊i_i, 0)
    b
end

"""
    pfPQ(; P=0, Q=0, name=:pqbus)

Create a PQ bus for power flow analysis.

A PQ bus has specified active and reactive power injections.
The voltage magnitude and phase angle are determined by the power flow solution.
"""
function pfPQ(; P=0, Q=0, name=:pqbus, kwargs...)
    @named pq = Library.PQConstraint(; P, Q)
    mtkbus = MTKBus(pq; name)
    b = compile_bus(mtkbus; kwargs...)
    set_voltage!(b; mag=1, arg=0)
    set_current!(b; P, Q)
    b
end

"""
    pfShunt(; B, G=0, name=:Shunt)

Create a bus with pure static shunt admittance for power flow analysis.
"""
function pfShunt(; B, G=0, name=:Shunt, kwargs...)
    @named shunt = Library.StaticShunt(; B, G)
    mtkbus = MTKBus(shunt; name)
    b = compile_bus(mtkbus; kwargs...)
    set_voltage!(b; mag=1, arg=0)
    set_current!(b, -(G + im*B))
    b
end

"""
    ispfmodel(cf::NetworkDynamics.ComponentModel)

Check if a component model is suitable for power flow analysis.

A component model is considered a valid power flow model if it has no dynamics,
i.e., either no states or a zero mass matrix.

## Returns
- `true` if the component is suitable for power flow analysis
- `false` otherwise
"""
function ispfmodel(cf::NetworkDynamics.ComponentModel)
    # no states or mass matrix 0
    if NetworkDynamics.dim(cf) == 0 || cf.mass_matrix==LinearAlgebra.UniformScaling(0)
        return true
    end
    return false
end

"""
    powerflow_model(cf::NetworkDynamics.ComponentModel; check=:error)

Extract or create a power flow component model from a dynamic component model.

1. If the component has `:pfmodel` metadata, use that model (after validation)
2. If the component is already a valid power flow model (i.e. no ODE, just constraints), return it as-is

## Parameters
- `cf`: The component model to extract/create a power flow model from
- `check`: Validation behavior (`:error`, `:warn`, or `:none`) when model validation fails

## Returns
- A component model suitable for power flow analysis (no dynamics)

## Validation
The returned model must satisfy [`ispfmodel`](@ref) criteria:
- Either no states or zero mass matrix (no dynamics)

See also: [`ispfmodel`](@ref), [`pfSlack`](@ref), [`pfPV`](@ref), [`pfPQ`](@ref)
"""
function powerflow_model(cf::NetworkDynamics.ComponentModel; check=:error)
    if has_pfmodel(cf)
        pfm = get_pfmodel(cf)
        if !ispfmodel(pfm)
            str = "Provided :pfmodel for :$(cf.name) is no valid powerflow model! You may overwrite this check by passing `check=:warn/:none`."
            check == :error && error(str)
            check == :warn && @warn str
        end
        return pfm
    end
    if !ispfmodel(cf)
        str = "Cannot create PF component model from :$(cf.name)! Please provide :pfmodel metadata! You may overwrite this check by passing `check=:warn/:none`."
        check == :error && error(str)
        check == :warn && @warn str
    end
    return cf
end

"""
    powerflow_model(nw::Network; check=:error)

Create a power flow network model from a dynamic network model.

This method applies [`powerflow_model`](@ref) to all vertex and edge components
in the network, creating a new network suitable for steady-state power flow analysis.

## Parameters
- `nw`: The network to create a power flow model from
- `check`: Validation behavior (`:error`, `:warn`, or `:none`) passed to component-level `powerflow_model` calls

## Returns
- A new `Network` with the same graph structure but power flow component models

See also: [`solve_powerflow`](@ref)
"""
function powerflow_model(nw::Network; check=:error)
    g = nw.im.g
    vfs = powerflow_model.(nw.im.vertexm; check);
    efs = powerflow_model.(nw.im.edgem; check);
    # graph elements like :bus1 => :bus2 might be wrong since the powerflow models
    # might have different names
    Network(g, vfs, efs; check_graphelement=false)
end

####
#### Power Flow Model Accessors
####

"""
    has_pfmodel(c::ComponentModel)
    has_pfmodel(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has a power flow model in metadata.

See also: [`get_pfmodel`](@ref), [`set_pfmodel!`](@ref).
"""
has_pfmodel(c::NetworkDynamics.ComponentModel) = has_metadata(c, :pfmodel)
has_pfmodel(nw::Network, idx::NetworkDynamics.VCIndex) = has_pfmodel(getcomp(nw, idx))
has_pfmodel(nw::Network, idx::NetworkDynamics.ECIndex) = has_pfmodel(getcomp(nw, idx))

"""
    get_pfmodel(c::NetworkDynamics.ComponentModel)
    get_pfmodel(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the power flow model for the component model.
May error if no power flow model is present. Use `has_pfmodel` to check first.

See also: [`has_pfmodel`](@ref), [`set_pfmodel!`](@ref).
"""
get_pfmodel(c::NetworkDynamics.ComponentModel) = get_metadata(c, :pfmodel)
get_pfmodel(nw::Network, idx::NetworkDynamics.VCIndex) = get_pfmodel(getcomp(nw, idx))
get_pfmodel(nw::Network, idx::NetworkDynamics.ECIndex) = get_pfmodel(getcomp(nw, idx))

"""
    set_pfmodel!(c::NetworkDynamics.ComponentModel, model)
    set_pfmodel!(nw::Network, idx::Union{VIndex,EIndex}, model)

Sets the power flow model for the component.
Overwrites any existing power flow model.

See also [`delete_pfmodel!`](@ref), [`get_pfmodel`](@ref).
"""
function set_pfmodel!(c::NetworkDynamics.ComponentModel, model)
    check_pfmodel(model, c)
    set_metadata!(c, :pfmodel, model)
end
set_pfmodel!(nw::Network, idx::NetworkDynamics.VCIndex, model) = set_pfmodel!(getcomp(nw, idx), model)
set_pfmodel!(nw::Network, idx::NetworkDynamics.ECIndex, model) = set_pfmodel!(getcomp(nw, idx), model)

"""
    delete_pfmodel!(c::NetworkDynamics.ComponentModel)
    delete_pfmodel!(nw::Network, idx::Union{VIndex,EIndex})

Removes the power flow model from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the power flow model existed and was removed, `false` otherwise.

See also: [`set_pfmodel!`](@ref).
"""
delete_pfmodel!(c::NetworkDynamics.ComponentModel) = delete_metadata!(c, :pfmodel)
delete_pfmodel!(nw::Network, idx::NetworkDynamics.VCIndex) = delete_pfmodel!(getcomp(nw, idx))
delete_pfmodel!(nw::Network, idx::NetworkDynamics.ECIndex) = delete_pfmodel!(getcomp(nw, idx))

"""
    solve_powerflow(nw::Network;
                    pfnw = powerflow_model(nw),
                    pfs0 = NWState(pfnw),
                    fill_busbar_defaults=true,
                    use_guesses=true,
                    sparse=nv(pfnw) > 50,
                    verbose=true,
                    kwargs...)

Solve the power flow equations for a given network.

Uses [`find_fixpoint`](@extref NetworkDynamics.find_fixpoint) from NetworkDynamics to solve the algebraic power flow equations.

## Parameters
- `nw`: The dynamic network model
- `pfnw`: The power flow network model (default: created from `nw`)
- `pfs0`: Initial state for the power flow calculation
- `fill_busbar_defaults`: Whether to fill missing default values for busbar states (i.e. u_r=1 u_i=0)
- `use_guesses`: Whether to fall back to "guess" values instead of "default" values if available.
- `sparse`: Whether to use a sparse solver (default: for networks with more than 50 buses)
- `verbose`: Whether to print the power flow solution
- `kwargs`: Additional keyword arguments passed to `find_fixpoint`

## Returns
- A `NWState` containing the solved power flow solution

See also [`initialize_from_pf`](@ref).
"""
function solve_powerflow(
    nw::Union{Network,Nothing};
    pfnw = powerflow_model(nw),
    pfs0 = NWState(pfnw),
    fill_busbar_defaults=true,
    use_guesses=true,
    sparse = nv(pfnw) > 50,
    verbose=true,
    alg = nothing,
    kwargs...
)
    # don't enforce this, check happens in `powerflow_model`
    # pfnw.mass_matrix == LinearAlgebra.UniformScaling(0) || error("Powerflow model must have a mass matrix of 0!")

    if isnothing(alg)
        if sparse && isnothing(pfnw.jac_prototype)
            alg = try
                set_jac_prototype!(pfnw)
                NonlinearSolve.FastShortcutNonlinearPolyalg(linsolve=NonlinearSolve.LinearSolve.KLUFactorization())
            catch e
                @warn "Could not set sparse jacobian prototype for powerflow model! Falling back to dense solver. Error: $e"
                NonlinearSolve.FastShortcutNonlinearPolyalg()
            end
        else
            alg = NonlinearSolve.FastShortcutNonlinearPolyalg()
        end
    end

    if fill_busbar_defaults && any(isnan, uflat(pfs0))
        urinds = generate_indices(pfnw, VIndex(:), :busbar₊u_r, s=true, obs=false, out=false, in=false, p=false)
        nanidx = findall(isnan, pfs0[urinds])
        if !isempty(nanidx)
            pfs0[urinds[nanidx]] .= 1.0
        end
        uiinds = generate_indices(pfnw, VIndex(:), :busbar₊u_i, s=true, obs=false, out=false, in=false, p=false)
        nanidx = findall(isnan, pfs0[uiinds])
        if !isempty(nanidx)
            pfs0[uiinds[nanidx]] .= 0.0
        end
    end
    if use_guesses && any(isnan, uflat(pfs0))
        for (i, idx) in enumerate(SII.variable_symbols(pfs0))
            isnan(uflat(pfs0)[i]) || continue
            has_guess(pfnw, idx) || continue
            uflat(pfs0)[i] = get_guess(pfnw, idx)
        end
    end

    pfs = find_fixpoint(pfnw, pfs0; alg, kwargs...)
    verbose && show_powerflow(pfs)

    return pfs
end

initialize_from_pf_docstring = raw"""
    initialize_from_pf[!](
        nw::Network;
        default_overrides=nothing,
        verbose = true,
        subverbose = false,
        pfnw = powerflow_model(nw),
        pfs0 = NWState(pfnw),
        pfs = solve_powerflow(nw; pfnw, pfs0, verbose),
        sparsepf = nv(nw) > 50,
        check = :error,
        kwargs...
    )

Initialize a dynamic network model from a power flow solution.

This function performs a two-step initialization process:
1. Solve the power flow equations for the network
2. Use the power flow solution to initialize the dynamic model

Per default, the powerflow solution is computed as

```julia
pfnw = powerflow_model(nw)                     # get powerflow network model
pfs0 = NWState(pfnw)                           # initial condition for network model
pfs = solve_powerflow(nw; pfnw, pfs0, verbose) # solve powerflow
```

You can override any of these steps by providing `pfnw`, `pfs0`, or `pfs` directly as a keyword argument.

There are two versions of this function: a mutating one (!-at the end of name) and a non-mutating version.
The mutating version uses `initialize_componentwise!` internally, the non-mutating one `initialize_componentwise`.
When the mutating version is used, `NWState(nw)` after initialization will return the same initialized
state again, as it is stored in the metadata.

## Parameters
- `nw`: The dynamic network model to initialize
- `default_overrides`: Is added *in addition* to the interface values extracted from the power flow solution. This allows you to change parameters for example.
- `verbose`: Whether to print information about the power flow solution (default: true)
- `subverbose`: Whether to print detailed information during component initialization (default: false). Can be Vector [VIndex(1), EIndex(3), ...] for selective output
- `pfnw`: Power flow network model (default: created from `nw` using `powerflow_model`)
- `pfs0`: Initial state for power flow calculation (default: created from `pfnw`)
- `pfs`: Power flow solution (default: calculated using `solve_powerflow`)
- `sparsepf`: Whether to use a sparse solver for power flow (default: for networks with more than 50 buses)
- `check`: Check per-unit base consistency by applying [`check_base_consistency`](@ref) to
  the initialized state. Either `:error` (default), `:warn` or `:none`.
- Additional keyword arguments are passed to `initialize_componentwise[!]`

## Returns
- A fully initialized network state

See also: [`solve_powerflow`](@ref), [`initialize_componentwise`](@extref NetworkDynamics.initialize_componentwise), [`interface_values`](@extref NetworkDynamics.interface_values)
"""
@doc initialize_from_pf_docstring
initialize_from_pf(nw; kw...) = _init_from_pf(initialize_componentwise, nw; kw...)
@doc initialize_from_pf_docstring
initialize_from_pf!(nw; kw...) = _init_from_pf(initialize_componentwise!, nw; kw...)
function _init_from_pf(
    initf, nw;
    verbose = true,
    subverbose = false,
    pfnw = nothing,
    pfs0 = nothing,
    pfs = nothing,
    sparsepf = nv(nw) > 50,
    check = :error,
    kwargs...
)
    if isnothing(pfs)
        pfnw = isnothing(pfnw) ? powerflow_model(nw) : pfnw
        pfs0 = isnothing(pfs0) ? NWState(pfnw) : pfs0
        pfs = solve_powerflow(nw; pfnw, pfs0, verbose, sparse=sparsepf)
    end

    interface_vals = interface_values(pfs)
    pfinitconstraints = specialize_pfinitconstraints(nw, pfs)
    pfinitformulas = specialize_pfinitformulas(nw, pfs)

    if haskey(kwargs, :additional_initconstraint)
        throw(ArgumentError("Cannot pass `additional_initconstraint` keyword down to \
            `initialize_componentwise[!]`! It is used internally to pass powerflow initconstraints. \
            If you need this functionality, please use `initialize_componentwise[!]` directly!"))
    end
    if haskey(kwargs, :additional_initformula)
        throw(ArgumentError("Cannot pass `additional_initformula` keyword down to \
            `initialize_componentwise[!]`! It is used internally to pass powerflow initformula. \
            If you need this functionality, please use `initialize_componentwise[!]` directly!"))
    end
    if haskey(kwargs, :default_overrides)
        default_overrides = merge(interface_vals, kwargs[:default_overrides])
        kwargs = filter(p -> p.first != :default_overrides, kwargs)
    else
        default_overrides = interface_vals
    end

    s0 = initf(
        nw;
        default_overrides,
        additional_initconstraint = pfinitconstraints,
        additional_initformula = pfinitformulas,
        verbose, subverbose, kwargs...
    )

    check_base_consistency(s0; check)

    s0
end

"""
    show_powerflow(s::NWState/Network)

Display power flow results in a tabular format.

Extract and format power flow solution data from a network state, showing bus-level information
including voltage magnitudes, phase angles, active power, and reactive power.
"""
function show_powerflow(s::NWState)
    NV = nv(extract_nw(s))
    dict = OrderedDict()
    dict["N"] = 1:NV
    dict["Bus Names"] = [cf.name for cf in extract_nw(s).im.vertexm]
    dict["vm [pu]"] = s[VIndex(1:NV, :busbar₊u_mag)]
    dict["varg [deg]"] = rad2deg.(s[VIndex(1:NV, :busbar₊u_arg)])
    dict["P [pu]"] = s[VIndex(1:NV, :busbar₊P)]
    dict["Q [pu]"] = s[VIndex(1:NV, :busbar₊Q)]

    DataFrame(dict)
end
show_powerflow(nw::Network) = show_powerflow(NWState(nw))

function check_pfmodel(pfm, m)
    in = NetworkDynamics.insym_flat(m)
    pfin = NetworkDynamics.insym_flat(pfm)
    out = NetworkDynamics.outsym_flat(m)
    pfout = NetworkDynamics.outsym_flat(pfm)
    if in != pfin || out != pfout
        @warn "Power flow model $(pfm.name) does not have the same input/output \
               variables as the original model $(m.name)! \n\
               - Model inputs:   $in\n\
               - PF mod inputs:  $pfin\n\
               - Model outputs:  $out\n\
               - PF mod outputs: $pfout\n\
               Maybe you're trying to attach an voltage source pf model to a current source component?"
    end
end

####
#### Per-unit base consistency
####
"""
    check_base_consistency(s::NWState; check=:error)

Verify that the per-unit bases of an initialized network agree with one another. Called
automatically at the end of [`initialize_from_pf`](@ref); `check` picks the severity
(`:error`, `:warn` or `:none`).

The *global* bases (Sbase, ωbase, ωframe) must be identical across every bus and every line.

The voltage base is genuinely per bus, so it is checked against the neighbors instead:

- each line end's `Vbase` must equal `busbar₊Vbase` of the bus it is attached to. This is
  per end, so a transformer with a different base on each side is fine — each end is
  compared against its own bus,
- each satellite (injector) bus must carry the same `busbar₊Vbase` as its hub, since it
  sits electrically *at* the hub bus.

Components which do not carry these symbols at all — a hand-rolled, purely per-unit model
without a [`SystemBase`](@ref), say — are skipped rather than failed; mixing based and
baseless components is legitimate. Missing and `NaN` values are skipped for the same reason
they are harmless: they fail loudly further downstream, rather than quietly producing a
plausible wrong number, which is the only failure mode this check exists to catch.
"""
function check_base_consistency(s::NWState; check=:error)
    check === :none && return nothing
    if check !== :error && check !== :warn
        throw(ArgumentError("`check` must be one of :error, :warn or :none, got $(repr(check))."))
    end
    nw = extract_nw(s)
    p = pflat(s)

    _vidxs(sym) = nv(nw) == 0 ? [] : SII.parameter_index(nw, VIndex(:, sym))
    _eidxs(sym) = ne(nw) == 0 ? [] : SII.parameter_index(nw, EIndex(:, sym))
    _val(i) = isnothing(i) || isnan(p[i]) ? nothing : p[i]

    msgs = String[]

    # global bases: exactly one value for the whole network, buses and lines alike
    allidxs = vcat(VIndex.(1:nv(nw)), EIndex.(1:ne(nw)))
    _globalvals(sym) = vcat(_val.(_vidxs(sym)), _val.(_eidxs(sym)))

    _check_uniform_base!(msgs, nw, :systembase₊Sbase,
        "power is mis-converted at every device/bus boundary",
        allidxs, _globalvals(:systembase₊Sbase))
    _check_uniform_base!(msgs, nw, :systembase₊ωbase,
        "time constants are silently rescaled per component",
        allidxs, _globalvals(:systembase₊ωbase))
    _check_uniform_base!(msgs, nw, :systembase₊ωframe,
        "the dq frames rotate apart, so no steady state exists",
        allidxs, _globalvals(:systembase₊ωframe))

    # `Vbase` is per bus, so it is checked locally. A single pass over the edges covers both
    # cases: a loopback edge attaches a satellite (src) to its hub (dst), any other edge
    # attaches two line ends to their respective buses.
    busV = _val.(_vidxs(:busbar₊Vbase))
    srcV = _val.(_eidxs(:src₊Vbase))
    dstV = _val.(_eidxs(:dst₊Vbase))
    for (e, edge) in pairs(nw.im.edgevec)
        if is_loopback(nw.im.edgem[e])
            _check_matching_base!(msgs, busV[edge.src], busV[edge.dst],
                "`busbar₊Vbase` of satellite $(_cdescr(nw, VIndex(edge.src)))",
                "its hub $(_cdescr(nw, VIndex(edge.dst)))",
                "the satellite sits electrically at the hub bus")
        else
            _check_matching_base!(msgs, srcV[e], busV[edge.src],
                "`src₊Vbase` of $(_cdescr(nw, EIndex(e)))",
                "its bus $(_cdescr(nw, VIndex(edge.src)))",
                "implies an unintended transformer ratio and wrong SI observables")
            _check_matching_base!(msgs, dstV[e], busV[edge.dst],
                "`dst₊Vbase` of $(_cdescr(nw, EIndex(e)))",
                "its bus $(_cdescr(nw, VIndex(edge.dst)))",
                "implies an unintended transformer ratio and wrong SI observables")
        end
    end

    isempty(msgs) && return nothing
    str = "The per-unit bases of this network are inconsistent:\n" * join(msgs, "\n") * "\n\
           Bases are copied into the components at *construction* time: set the global ones \
           with `set_Sbase!`/`set_ωbase!` before building, or explicitly via `SystemBase(…)`, \
           and the per-bus voltage base via `MTKBus(…; Vbase)`/`compile_bus(…; Vbase)`. \
           Pass `check=:warn` or `check=:none` to `initialize_from_pf` to bypass this check."
    check === :error ? error(str) : @warn str
    nothing
end
function _check_uniform_base!(msgs, nw, sym, consequence, idxs, vals)
    keep = .!isnothing.(vals)
    any(keep) || return msgs
    idxs = idxs[keep]
    vals = convert(Vector{Float64}, vals[keep])
    allequal(vals) && return msgs

    groups = [(v, idxs[vals .== v]) for v in unique(vals)]
    sort!(groups; by=g -> length(g[2]), rev=true)
    parts = map(enumerate(groups)) do (i, (v, gidxs))
        if i == 1 && length(gidxs) > 1
            "$v on $(length(gidxs)) components"
        else
            "$v on $(_cdescrs(nw, gidxs))"
        end
    end
    push!(msgs, " - `$sym` differs across the network ($consequence): " * join(parts, ", "))
end
function _check_matching_base!(msgs, a, b, adescr, bdescr, consequence)
    (isnothing(a) || isnothing(b) || a == b) && return msgs
    push!(msgs, " - $adescr is $a but $bdescr has $b ($consequence)")
end
_cdescr(nw, idx) = "$idx ($(nw[idx].name))"
function _cdescrs(nw, idxs; max=8)
    shown = join((_cdescr(nw, idx) for idx in Iterators.take(idxs, max)), ", ")
    length(idxs) > max ? shown * " and $(length(idxs) - max) more" : shown
end
