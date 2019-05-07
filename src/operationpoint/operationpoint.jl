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
