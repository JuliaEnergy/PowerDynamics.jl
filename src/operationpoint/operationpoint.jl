using NLsolve: nlsolve, converged

struct RootRhs
    rhs
end
function (rr::RootRhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    dx
end

function RootRhs(of::ODEFunction)
    RootRhs(of.f)
end

function find_operationpoint(pg::PowerGrid, ic_guess = nothing)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        throw(OperationPointError("currently not making any checks concerning assumptions of whether its possible to find the fixed point"))
    end
    if SwingEq ∈ pg.nodes .|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end

    if ic_guess === nothing
        system_size = systemsize(pg)
        ic_guess = ones(system_size)
    end

    rr = RootRhs(ode_function(pg))
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end
