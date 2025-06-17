function pfSlack(; V=missing, δ=missing, u_r=missing, u_i=missing)
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

    mtkbus = MTKBus(slack, name=:slackbus)
    b = Bus(mtkbus)
    set_voltage!(b, u_r + im * u_i)
    b
end

function pfPV(; P, V)
    @named pv = Library.PVConstraint(; P, V)
    mtkbus = MTKBus(pv, name=:pvbus)
    b = Bus(mtkbus)
    set_voltage!(b; mag=V, arg=0)
    b
end

function pfPQ(; P=0, Q=0)
    @named pq = Library.PQConstraint(; P, Q)
    mtkbus = MTKBus(pq, name=:pqbus)
    b = Bus(mtkbus)
    set_voltage!(b; mag=1, arg=0)
    b
end

function ispfmodel(cf::NetworkDynamics.ComponentModel)
    # cannot have freep nor freeu
    isempty(freep(cf)) || return false
    isempty(freeu(cf)) || return false
    # no states or mass matrix 0
    if NetworkDynamics.dim(cf) == 0 || cf.mass_matrix==LinearAlgebra.UniformScaling(0)
        return true
    end
    return false
end

function powerflow_model(cf::NetworkDynamics.ComponentModel)
    if has_metadata(cf, :pfmodel)
        pfm = get_metadata(cf, :pfmodel)
        if !ispfmodel(pfm)
            error("Provided :pfmodel for :$(cf.name) is no valid powerflow model!")
        end
        return pfm
    elseif ispfmodel(cf)
        return cf
    elseif cf isa VertexModel && has_default(cf, :busbar₊u_r) && has_default(cf, :busbar₊u_i)
        @warn "No powerflow model given for :$(cf.name), using slack with default voltage!"
        return pfSlack(u_r=get_default(cf, :busbar₊u_r), u_i=get_default(cf, :busbar₊u_i))
    end
    error("Cannot create PF component model from :$(cf.name)! Please proved :pfmodel metadata!")
end

function powerflow_model(nw::Network)
    g = nw.im.g
    vfs = powerflow_model.(nw.im.vertexm);
    efs = powerflow_model.(nw.im.edgem);
    Network(g, vfs, efs)
end

"""
    solve_powerflow(nw::Network;
                    pfnw = powerflow_model(nw),
                    pfs0 = NWState(nw),
                    verbose=true)

Solve the power flow equations for a given network.

Uses [`find_fixpoint`](@extref) from NetworkDynamics to solve the algebraic power flow equations.

## Parameters
- `nw`: The dynamic network model
- `pfnw`: The power flow network model (default: created from `nw`)
- `pfs0`: Initial state for the power flow calculation
- `verbose`: Whether to print the power flow solution

## Returns
- A `NWState` containing the solved power flow solution

See also [`intitialize_from_pf`](@ref).
"""
function solve_powerflow(
    nw::Network;
    pfnw = powerflow_model(nw),
    pfs0 = NWState(nw),
    verbose=true
)
    pfnw.mass_matrix == LinearAlgebra.UniformScaling(0) || error("Powerflow model must have a mass matrix of 0!")

    pfs0 = NWState(pfnw)
    uf = uflat(pfs0)
    pf = pflat(pfs0)
    any(isnan, uf) && error("Initial state for powerflow model contains NaNs!")
    any(isnan, pf) && error("Parameters for powerflow model contain NaNs!")

    pfs = find_fixpoint(pfnw, pfs0)
    verbose && show_powerflow(pfs)

    return pfs
end

initialize_from_pf_docstring = raw"""
    initialize_from_pf[!](
        nw::Network;
        verbose = true,
        subverbose = false,
        pfnw = powerflow_model(nw),
        pfs0 = NWState(pfnwnw),
        pfs = solve_powerflow(pfnw, pfs0; verbose=verbose),
        kwargs...
    )

Initialize a dynamic network model from a power flow solution.

This function performs a two-step initialization process:
1. Solve the power flow equations for the network
2. Use the power flow solution to initialize the dynamic model

There are two versions of this function: a mutating one (!-at the end of name) and a non-mutating version.
The mutating version uses `initialize_componentwise!` internally, the non-mutating one `initialize_componentwise`.
When the mutating version is used, `NWState(nw)` after initialization will return the same initialized
state again, as it is stored in the metadata.

## Parameters
- `nw`: The dynamic network model to initialize
- `verbose`: Whether to print information about the power flow solution (default: true)
- `subverbose`: Whether to print detailed information during component initialization (default: false)
- `pfnw`: Power flow network model (default: created from `nw` using `powerflow_model`)
- `pfs0`: Initial state for power flow calculation (default: created from `pfnw`)
- `pfs`: Power flow solution (default: calculated using `solve_powerflow`)
- Additional keyword arguments are passed to `initialize_componentwise[!]`

## Returns
- A fully initialized network state

See also: [`solve_powerflow`](@ref), [`initialize_componentwise`](@extref), [`interface_values`](@extref)
"""
@doc initialize_from_pf_docstring
initialize_from_pf(nw; kw...) = _init_from_pf(initialize_componentwise, nw; kw...)
@doc initialize_from_pf_docstring
initialize_from_pf!(nw; kw...) = _init_from_pf(initialize_componentwise!, nw; kw...)
function _init_from_pf(
    initf, nw;
    verbose = true,
    subverbose = false,
    pfnw = powerflow_model(nw),
    pfs0 = NWState(pfnwnw),
    pfs = solve_powerflow(pfnw, pfs0; verbose=verbose),
    kwargs...
)
    interface_vals = interface_values(pfs)
    initf(nw; verbose, subverbose, kwargs...)
end

show_powerflow(nw::Network) = show_powerflow(NWState(nw))
function show_powerflow(s::NWState)
    NV = nv(extract_nw(s))
    dict = OrderedDict()
    dict["N"] = 1:NV
    dict["Bus Names"] = [cf.name for cf in extract_nw(s).im.vertexm]
    dict["vm [pu]"] = s[vidxs(1:NV, :busbar₊u_mag)]
    dict["varg [deg]"] = rad2deg.(s[vidxs(1:NV, :busbar₊u_arg)])
    dict["P [pu]"] = s[vidxs(1:NV, :busbar₊P)]
    dict["Q [pu]"] = s[vidxs(1:NV, :busbar₊Q)]

    DataFrame(dict)
end
