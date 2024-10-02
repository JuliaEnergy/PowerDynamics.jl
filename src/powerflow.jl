function pfSlack(; V, δ=0)
    @named slack = Library.VδConstraint(; V, δ)
    mtkbus = MTKBus(slack, name=:slackbus)
    b = Bus(mtkbus)
    set_voltage!(b; mag=V, arg=δ)
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

function ispfmodel(cf::NetworkDynamics.ComponentFunction)
    isempty(freep(cf)) || return false
    NetworkDynamics.isstatic(cf) && return true
    isempty(freeu(cf)) || return false
    cf.mass_matrix==LinearAlgebra.UniformScaling(0) && return true
    return false
end

function powerflow_model(cf::NetworkDynamics.ComponentFunction)
    if NetworkDynamics.has_metadata(cf, :pfmodel)
        pfm = NetworkDynamics.get_metadata(cf, :pfmodel)
        if !ispfmodel(pfm)
            error("Provided :pfmodel for $(cf.name) is no valid powerflow model!")
        end
        return pfm
    elseif ispfmodel(cf)
        return cf
    end
    error("Cannot create PF component model from $(cf.name)! Please proved :pfmodel metadata!")
end

function powerflow_model(nw::Network)
    g = nw.im.g
    vfs = powerflow_model.(nw.im.vertexf);
    efs = powerflow_model.(nw.im.edgef);
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
    # TODO: get input states as observables
    _, aggbuf = NetworkDynamics.get_ustacked_buf(pfs)

    for i in 1:nv(nw)
        set_voltage!(nw.im.vertexf[i], pfs.v[i, :busbar₊u_r] + im * pfs.v[i, :busbar₊u_i])
        ir, ii = aggbuf[nw.im.v_aggr[i]]
        set_current!(nw.im.vertexf[i], ir + im*ii)
    end

    println("Found powerflow solution:")
    df = show_powerflow(nw)
    println(df)
    nothing
end

function show_powerflow(nw::Network)
    # df = DataFrame()
    dict = OrderedDict()
    dict["N"] = 1:nv(nw)
    dict["Bus Names"] = [cf.name for cf in nw.im.vertexf]
    s = NWState(nw)
    u = Complex.(s.v[:, :busbar₊u_r], s.v[:, :busbar₊u_i])
    _, aggbuf = NetworkDynamics.get_ustacked_buf(s)
    i = [Complex(aggbuf[nw.im.v_aggr[k]]...) for k in 1:nv(nw)]
    S = u .* conj.(i)

    dict["vm [pu]"] = abs.(u)
    dict["varg [deg]"] = rad2deg.(angle.(u))
    dict["P [pu]"] = real.(S)
    dict["Q [pu]"] = imag.(S)

    DataFrame(dict)
end

function initialize!(nw::Network; verbose=true)
    for cf in nw.im.vertexf
        NetworkDynamics.initialize_component!(cf; verbose=false)
        res = LinearAlgebra.norm(get_metadata(cf, :init_residual))
        if verbose
            if res < 1e-8
                printstyled("$(cf.name) successful! (residual=$res)\n", color=:green)
            else
                printstyled("$(cf.name) failed! (residual=$res)\n", color=:red)
            end
        else
            res > 1e-8 && @warn "Initialization of $(cf.name) failed! Residual: $res"
        end

    end
    nothing
end
