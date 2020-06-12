using NLsolve: nlsolve, converged
using SteadyStateDiffEq
using OrdinaryDiffEq: ODEProblem, Rodas5

"""
```Julia
initial_guess(pg)
```
makes a type-specific initial guess to help the operation point search

# Arguments
- `pg`: a [`PowerGrid`](@ref) instance


The voltage of all nodes is guessed as the slack voltage in the system.
See also the documentation of [`guess`](@ref).
"""
function initial_guess(pg)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers. Default voltage guess: u = 1 + 0j [pu]."
        voltage_guess = complex(1.0, 0.0)
    else
        sl = findfirst(SlackAlgebraic ∈ pg.nodes .|> typeof)
        slack = pg.nodes[sl]
        voltage_guess = slack.U
    end

    return initial_guess(pg, voltage_guess)
end

"""
```Julia
initial_guess(pg, complex_voltages)
```
makes a type-specific initial guess to help the operation point search

# Arguments
- `pg`: a [`PowerGrid`](@ref) instance
- `complex_voltages`: voltage guess, either scalar or vector of `length(pg.nodes)`


The voltage of all nodes is guessed as the respective entry in `complex_voltages`,
the latter can e.g. be the result of a power flow solution.
See also the documentation of [`guess`](@ref).
"""
function initial_guess(pg, complex_voltages)
    type_guesses = guess.(pg.nodes, complex_voltages)
    return vcat(type_guesses...)
end


"""
```Julia
guess(node, voltage_guess)
```

creates an initial guess for an operation point

# Arguments
- `node`: a node in a [`PowerGrid`](@ref), subtype of [`AbstractNode`](@ref)
- `voltage_guess`: guess for the node's voltage, defaults to slack set point

The guesses are specific to the node types. In the default case, the
node voltage set to the `voltage_guess` while internal variables are
initialised at 0.

It is possible to add methods for user-defined node types.

# Examples

+ guess(::Type{SlackAlgebraic}, U) = [real(U), imag(U)]         #[:u_r, :u_i]
+ guess(::Type{PQAlgebraic}, U) = [real(U), imag(U)]            #[:u_r, :u_i]
+ guess(::Type{ThirdOrderEq}, U) = [real(U), imag(U), 0., 0.]   #[:u_r, :u_i, :θ, :ω]
+ guess(::Type{<:AbstractNode}, U) = [real(U), imag(U), 0., ..., 0.]   #[:u_r, :u_i, ...]
"""
function guess end

# the default case
function guess(n::AbstractNode, voltage_guess)
    state = zeros(dimension(n))
    state[1] = real(voltage_guess)
    state[2] = imag(voltage_guess)
    return state
end

# In case the system has a slack bus, make sure the voltage_guess is compatible.
function guess(n::SlackAlgebraic, voltage_guess)
    @assert n.U ≈ voltage_guess
    return [real(n.U), imag(n.U)]
end

"""
```Julia
find_operationpoint(pg, ic_guess, p0, t0)
```
returns an operation point

# Arguments
- `pg`: a [`PowerGrid`](@ref) instance

# Optional Arguments
- `ic_guess`: initial guess for the operation point
- `p0`: possibility to pass parameters to `rhs(pg)`
- `t0`: specify evaluation time point for non-autonomus dynamics (defaults to t0 = 0)
"""
function find_operationpoint(
    pg::PowerGrid,
    ic_guess = nothing;
    p0 = nothing,
    t0 = 0.0,
    sol_method = :rootfind,
    solver_kwargs...
)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find a operation point."
    end
    if SwingEq ∈ pg.nodes .|> typeof
        throw(OperationPointError("Found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)."))
    end

    if ic_guess === nothing
        ic_guess = initial_guess(pg)
    end

    if sol_method == :nlsolve
        return _find_operationpoint_nlsolve(pg, ic_guess, p0, t0; solver_kwargs...)
    elseif sol_method == :rootfind
        return _find_operationpoint_rootfind(pg, ic_guess, p0, t0; solver_kwargs...)
    elseif sol_method == :steadystate
        return _find_operationpoint_steadystate(pg, ic_guess, p0, t0; solver_kwargs...)
    else
        throw(OperationPointError("$sol_method is not supported. Pass either `:nlsolve`, `:rootfind` or `:steadystate`"))
    end
end

function _find_operationpoint_steadystate(pg, ic_guess, p0, t0; kwargs...)
    ode = rhs(pg)
    op_prob = ODEProblem(ode, ic_guess, Inf)
    sol = solve(
        SteadyStateProblem(op_prob),
        DynamicSS(Rodas5(); kwargs...),
    )
    if sol.retcode == :Success
        return State(pg, sol.u)
    else
        throw(OperationPointError("The operation point search did not converge."))
    end
end

function _find_operationpoint_rootfind(pg, ic_guess, p0, t0; kwargs...) #solver=Rodas5(), abstol = 1e-8, reltol = 1e-6, tspan = Inf)
    ode = rhs(pg)
    op_prob = ODEProblem(ode, ic_guess, Inf)
    sol = solve(
        SteadyStateProblem(op_prob),
        SSRootfind(;kwargs...),
    )
    if sol.retcode == :Success
        return State(pg, sol.u)
    else
        throw(OperationPointError("The operation point search did not converge."))
    end
end

function _find_operationpoint_nlsolve(pg, ic_guess, p0, t0; kwargs...)
    # construct rhs for rootfinding
    rhs_pg = rhs(pg)
    rr = (dx, x) -> rhs_pg(dx, x, p0, t0)

    nl_res = nlsolve(rr, ic_guess; kwargs...)

    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        throw(OperationPointError("The operation point search did not converge."))
    end
end
