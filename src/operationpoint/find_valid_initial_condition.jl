using OrdinaryDiffEq: ODEFunction
using NLsolve: nlsolve, converged
using LinearAlgebra: pinv, UniformScaling

struct RootRhs_ic
    rhs
    mpm
end
function (rr::RootRhs_ic)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mpm * dx .- dx
end

function RootRhs_ic(of::ODEFunction)
    if of.mass_matrix isa UniformScaling
        n = length(of.syms)
        mm = Array(of.mass_matrix, n, n)
    else
        mm = of.mass_matrix
    end
    @assert mm !== nothing
    mpm = pinv(mm) * mm
    RootRhs_ic(of.f, mpm)
end


function find_valid_initial_condition(pg::PowerGrid, ic_guess)
    rr = RootRhs_ic(rhs(pg))
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return nl_res.zero #State(pg, nl_res.zero)
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end

export RootRhs_ic
export find_valid_initial_condition
