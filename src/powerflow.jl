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

function solve_powerflow!(nw::Network; verbose=true)
    pfnw = powerflow_model(nw)
    pfnw.mass_matrix == LinearAlgebra.UniformScaling(0) || error("Powerflow model must have a mass matrix of 0!")

    u0 = NWState(pfnw)
    uf = uflat(u0)
    pf = pflat(u0)
    any(isnan, uf) && error("Initial state contains NaNs!")
    any(isnan, pf) && error("Parameters contain NaNs!")
    prob = NonlinearProblem((du,u,p) -> pfnw(du,u,p,NaN), uf, pf)
    sol = solve(prob)
    if !SciMLBase.successful_retcode(sol.retcode)
        error("Powerflow did not converge! Retcode $(sol.retcode)")
    end

    pfs = NWState(pfnw, sol.u, pf)
    for i in 1:nv(nw)
        set_voltage!(nw.im.vertexm[i], pfs.v[i, :busbar₊u_r] + im * pfs.v[i, :busbar₊u_i])
        set_current!(nw.im.vertexm[i], pfs.v[i, :busbar₊i_r] + im * pfs.v[i, :busbar₊i_i])
    end

    show_powerflow(nw)
end

function show_powerflow(nw::Network)
    # df = DataFrame()
    dict = OrderedDict()
    dict["N"] = 1:nv(nw)
    dict["Bus Names"] = [cf.name for cf in nw.im.vertexm]
    # s = NWState(nw)
    u = Vector{Complex{Float64}}(undef, nv(nw))
    S = Vector{Complex{Float64}}(undef, nv(nw))

    try
        for (i, cf) in pairs(nw.im.vertexm)
            u[i] = get_voltage(cf)
            S[i] = get_power(cf)
        end
    catch e
        throw(ArgumentError("Could not extract voltage and power from vertex
            models, make sure that all bus models have default voltage and current
            set. For example by calling `solve_powerflow!`"))
    end

    dict["vm [pu]"] = abs.(u)
    dict["varg [deg]"] = rad2deg.(angle.(u))
    dict["P [pu]"] = real.(S)
    dict["Q [pu]"] = imag.(S)

    DataFrame(dict)
end

function initialize!(nw::Network; verbose=true)
    for cf in nw.im.vertexm
        fp = length(freep(cf))
        fu = length(freeu(cf))
        try
            NetworkDynamics.initialize_component!(cf; verbose=false)
        catch e
            println(e.msg)
            set_metadata!(cf, :init_residual, Inf)
        end
        res = LinearAlgebra.norm(get_metadata(cf, :init_residual))
        if verbose
            if res < 1e-8
                printstyled("$(cf.name) successful! ($fu/$fp free states/p; residual=$res)\n", color=:green)
            else
                printstyled("$(cf.name) failed! ($fu/$fp free states/p; residual=$res)\n", color=:red)
            end
        else
            res > 1e-8 && @warn "Initialization of $(cf.name) failed! Residual: $res"
        end

    end
    nothing
end
