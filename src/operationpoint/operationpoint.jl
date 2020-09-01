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
    if SlackAlgebraic ∉ collect(values(pg.nodes)) .|> typeof
        @warn "There is no slack bus in the system to balance powers. Default voltage guess: u = 1 + 0j [pu]."
        voltage_guess = complex(1.0, 0.0)
    else
        sl = findfirst(collect(values(pg.nodes).|> typeof).== SlackAlgebraic)
        bus_array=collect(values(pg.nodes))
        slack = bus_array[sl]
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
    bus_array=collect(values(pg.nodes))
    type_guesses = guess.(bus_array, complex_voltages)
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

# Argument
- `pg`: a [`PowerGrid`](@ref) instance

# Optional Arguments
- `ic_guess`: custom initial guess for the operation point

# Optional Keyword Arguments
- `p0`: possibility to pass parameters to `rhs(pg) (defaults to p0 = nothing)`
- `t0`: specify evaluation time point for non-autonomus dynamics (defaults to t0 = 0)
- `sol_method` : method used to locate an operation point, see comment below
- `sol_kwargs` : optional arguments passed to the solver

# Operation Point Solvers
- `:nlsolve` : uses `nlsolve` from `NLsolve` to find a root of the RHS with the trust region method (others can be selected)
- `:rootfind` : (default) uses `SSRootfind` from `SteadyStateDiffEq` to perform a NLsolve-based rootfind
- `:dynamic` : uses `DynamicSS` from `SteadyStateDiffEq` to integrate the dynamical equations until a steady state is reached

See also the documentation of the packages `SteadyStateDiffEq` and `NLsolve` to
see further finetuning options that can be passed via the `sol_kwargs` argument.

"""
function find_operationpoint(
    pg::PowerGrid,
    ic_guess = nothing;
    p0 = nothing,
    t0 = 0.0,
    sol_method = :rootfind,
    sol_kwargs...
)
    if SlackAlgebraic ∉ collect(values(pg.nodes)) .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find a operation point."
    end
    if SwingEq ∈ collect(values(pg.nodes)) .|> typeof
        throw(OperationPointError("Found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)."))
    end

    if ic_guess === nothing
        ic_guess = initial_guess(pg)
    end

    if sol_method == :nlsolve
        return _find_operationpoint_nlsolve(pg, ic_guess, p0, t0; sol_kwargs...)
    elseif sol_method == :rootfind
        return _find_operationpoint_rootfind(pg, ic_guess, p0, t0; sol_kwargs...)
    elseif sol_method == :dynamic
        return _find_operationpoint_steadystate(pg, ic_guess, p0, t0; sol_kwargs...)
    else
        throw(OperationPointError("$sol_method is not supported. Pass either `:nlsolve`, `:rootfind` or `:dynamic`"))
    end
end

function _find_operationpoint_steadystate(pg, ic_guess, p0, t0; kwargs...)
    ode = rhs(pg)
    op_prob = ODEProblem(ode, ic_guess, (t0, Inf), p0)
    sol = solve(
        SteadyStateProblem(op_prob),
        DynamicSS(Rodas5(); kwargs...),
    )
    if sol.retcode == :Success
        return State(pg, sol.u)
    else
        @warn "The operation point search did not converge. (dynamic method, $(sol.retcode))\n Fallback to rootfinding."
        _find_operationpoint_nlsolve(pg, ic_guess, p0, t0)
    end
end

function _find_operationpoint_rootfind(pg, ic_guess, p0, t0; kwargs...)
    ode = rhs(pg)
    op_prob = ODEProblem(ode, ic_guess, (t0, Inf), p0)
    sol = solve(
        SteadyStateProblem(op_prob),
        SSRootfind(;kwargs...),
    )
    if sol.retcode == :Success
        return State(pg, sol.u)
    else
        throw(OperationPointError("The operation point search did not converge. (rootfind method, $(sol.retcode))"))
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
        throw(OperationPointError("The operation point search did not converge. (nlsolve method)"))
    end
end
