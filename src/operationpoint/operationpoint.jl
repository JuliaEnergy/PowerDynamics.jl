using NLsolve: nlsolve, converged

struct RootRhs
    rhs
end
function (rr::RootRhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    dx
end

function find_operationpoint(pg::PowerGrid, ic_guess = nothing)
    if SlackAlgebraic ∉ collect(values(pg.nodes)) .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find the fixed point"
    end
    if SwingEq ∈ collect(values(pg.nodes)).|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end

    if ic_guess === nothing
        system_size = systemsize(pg)
        ic_guess = ones(system_size)
    end

    rr = RootRhs(rhs(pg))
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        throw(OperationPointError("Failed to find initial conditions on the constraint manifold!"))
    end
end
