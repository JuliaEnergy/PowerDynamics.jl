using NLsolve

struct RootRhs2
    rhs
end
function (rr::RootRhs2)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    dx
end

function RootRhs2(of::ODEFunction)
    RootRhs2(of.f)
end

function find_operationpoint(of::ODEFunction, ic_guess)
    rr = RootRhs2(of)
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end

struct constant_rotation_root
    nd # from NetworkDynamics with PowerDynamics conventions
    p
end

function (crr::constant_rotation_root)(x)
    all_v_idx = crr.nd.f.v_idx
    mm = crr.nd.mass_matrix
    v_idx = filter(idx -> mm[idx[1],idx[1]] > 0. &&  mm[idx[2],idx[2]] > 0., all_v_idx)

    dx = similar(x)
    crr.nd(dx, x, crr.p, 0.)

    dv_p = im*0.
    for idx in v_idx
        dv_p += complex(dx[idx[1]], dx[idx[2]]) / complex(x[idx[1]], x[idx[2]])
    end

    dv_p /= length(v_idx)

    i_omega = im * imag(dv_p)

    for idx in v_idx
        v = complex(x[idx[1]], x[idx[2]])
        res = complex(dx[idx[1]], dx[idx[2]]) - i_omega * v
        dx[idx[1]], dx[idx[2]] = real(res), imag(res)
    end


    # We want to set the quantity d(v_p^* v)/dt to zero rather than dv/dt at
    # dynamical nodes.
    # This is so that a homogeneous rotation of all voltage vectors together
    # is acceptable.

    #=
    v_p_conj = conj(complex(x[v_idx[end][1]], x[v_idx[end][2]]))
    dv_p_conj = conj(complex(dx[v_idx[end][1]], dx[v_idx[end][2]]))
    for idx in v_idx
        v = complex(x[idx[1]], x[idx[2]])
        dv = complex(dx[idx[1]], dx[idx[2]])
        res = dv_p_conj * v + v_p_conj * dv
        v_p_conj = conj(v)
        dv_p_conj = conj(dv)
        dx[idx[1]], dx[idx[2]] = real(res), imag(res)
    end
    =#
    #TODO Dynamic edges require further work here to take care of the phase
    # variation in the complex currents. It might be sensible to replace the
    # di/dt by the derivative of the power.
    dx
end

export find_operationpoint2

function find_operationpoint2(pg::PowerGrid, ic_guess)
    # TODO Check in pg.vertices for slack bus
    crr = constant_rotation_root(pg.network_dynamics, nothing)
    nl_res = nlsolve(crr, ic_guess)
    println("1")
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
        return nl_res.zero
    end
end
