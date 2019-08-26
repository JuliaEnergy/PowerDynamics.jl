using Pkg
Pkg.instantiate()
cd(@__DIR__)
using Plots
using NLsolve: nlsolve, converged
using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!,savevalues!

struct RootRhs
    rhs
end
function (rr::RootRhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    dx
end


# this doesn't work & hangs in an infinite loop
function simulate(rhs1, rhs2, x0, timespan,tspan_fault)

    problem = ODEProblem{true}(rhs1, x0, timespan)
    integrator = init(problem, Rodas4(autodiff=false))

    step!(integrator, tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs2

    step!(integrator, tspan_fault[2], true)

    # update integrator, clear error
    integrator.f = rhs1

    solve!(integrator)

    return integrator.sol
end

# this works, but leads to 3 seaprate solutions that need to be combined later
function simulate2(rhs1, rhs2, x0, timespan,tspan_fault)

    problem1 = ODEProblem{true}(rhs1,x0.vec,(timespan[1],tspan_fault))
    sol1 = solve(problem1, Rodas4(autodiff=false))

    final_state1 = sol1(:final)
    problem2 = ODEProblem{true}(rhs2,x0.vec,final_state1, (tspan_fault[1], tspan_fault[2]))
    sol2 = solve(problem2, Rodas4(autodiff=false))

    final_state2 = sol2(:final)
    # solve after clearance of the fault
    problem3 = ODEProblem{true}(rhs1,x0.vec,final_state2, (tspan_fault[2], last(timespan)))
    sol3 = solve(problem2, Rodas4(autodiff=false))

    return (sol1, sol2, sol3)
end

# this works, but the code looks hacky. Can't we do any better?
function simulate3(rhs1, rhs2, x0, timespan,tspan_fault)
    problem2 = ODEProblem{true}(rhs2,x0,(first(timespan),tspan_fault[2]))
    fault_integrator = init(problem2, Rodas4(autodiff=false))

    reinit!(fault_integrator, x0, t0=tspan_fault[1], tf=tspan_fault[2], erase_sol=false)
    savevalues!(fault_integrator)
    solve!(fault_integrator)

    problem1 = ODEProblem{true}(rhs1, fault_integrator.u, (tspan_fault[2], last(timespan)))
    integrator = init(problem1, Rodas4(autodiff=false))

    # Now the trick: copy solution object to new integrator and
    # make sure the counters are updated, otherwise sol is overwritten in the
    # next step.
    integrator.sol = fault_integrator.sol
    integrator.saveiter = fault_integrator.saveiter
    integrator.saveiter_dense = fault_integrator.saveiter_dense
    integrator.success_iter = fault_integrator.success_iter

    solve!(integrator)

    return integrator.sol
end


function find_operationpoint(rr::RootRhs, ic_guess = nothing)

    system_size = 6
    ic_guess = ones(system_size)

    nl_res = nlsolve(rr, ic_guess)
    print(nl_res)
    if converged(nl_res) == true
        return nl_res.zero
    else
        throw("Failed to find initial conditions on the constraint manifold!")
    end
end

begin
    P_1 =1
    H_1 =5
    P_2 = -1
    H_2 =5
    D_1=0.1
    D_2=1
    Ω_H_1 = 2*π*50/H_1
    Ω_H_2 = 2*π*50/H_2
    Γ=0.1
    V=1.
end

admittance = -1im/0.02 #+ 0.03
Y_matrix = [-admittance admittance;
            admittance -admittance]

function rhs!(dx, x, p, t)
    u_1 = x[1] + x[2] * im
    ω_1 = x[3]
    u_2 = x[4] + x[5] * im
    ω_2 = x[6]
    u = [u_1;u_2]
    i = Y_matrix*u
    i_1 = i[1]
    i_2 = i[2]
    p = real(u .* conj(i))
    dϕ_1 = ω_1
    dω_1 = ((P_1 - D_1 * ω_1) - p[1]) * Ω_H_1
    v_1 = abs(u_1)
    dv_1 = -Γ * (v_1 - V)
    du_1 = (u_1 / v_1) * dv_1 + u_1 * im * dϕ_1

    dϕ_2 = ω_2
    dω_2 = ((P_2 - D_2 * ω_2) - p[2]) * Ω_H_2
    v_2 = abs(u_2)
    dv_2 = -Γ * (v_2 - V)
    du_2 = (u_2 / v_2) * dv_2 + u_2 * im * dϕ_2
    try
        dx[1] = real(du_1)
        dx[2] = imag(du_1)
        dx[3] = dω_1
        dx[4] = real(du_2)
        dx[5] = imag(du_2)
        dx[6] = dω_2
        return nothing
    catch e
           if typeof(e) === UndefVarError
               throw(NodeDynamicsError("you need to provide $(e.var)"))
           else
               throw(e)
           end
       end
    end
P_1=0.9
function rhs2!(dx, x, p, t)
    u_1 = x[1] + x[2] * im
    ω_1 = x[3]
    u_2 = x[4] + x[5] * im
    ω_2 = x[6]
    u = [u_1;u_2]
    i = Y_matrix*u
    i_1 = i[1]
    i_2 = i[2]
    p = real(u .* conj(i))
    dϕ_1 = ω_1
    dω_1 = ((P_1 - D_1 * ω_1) - p[1]) * Ω_H_1
    v_1 = abs(u_1)
    dv_1 = -Γ * (v_1 - V)
    du_1 = (u_1 / v_1) * dv_1 + u_1 * im * dϕ_1

    dϕ_2 = ω_2
    dω_2 = ((P_2 - D_2 * ω_2) - p[2]) * Ω_H_2
    v_2 = abs(u_2)
    dv_2 = -Γ * (v_2 - V)
    du_2 = (u_2 / v_2) * dv_2 + u_2 * im * dϕ_2
    try
        dx[1] = real(du_1)
        dx[2] = imag(du_1)
        dx[3] = dω_1
        dx[4] = real(du_2)
        dx[5] = imag(du_2)
        dx[6] = dω_2
        return nothing
    catch e
           if typeof(e) === UndefVarError
               throw(NodeDynamicsError("you need to provide $(e.var)"))
           else
               throw(e)
           end
       end
    end

rr1=RootRhs(rhs!)
#rr2 = RootRhs(rhs2!)
operationpoint = find_operationpoint(rr1)

tspan_fault = (0.1,1.)
disturbed_node=2
result = simulate(rhs!,rhs2!,
    operationpoint,
    tspan_fault,
    (0., 1.))

#result = simulate(Perturbation(disturbed_node, :ω, Inc(0.5)), powergrid, operationpoint, timespan = (0.0,1.))
#result = simulate(LineFault(2,3),powergrid,operationpoint,timespan=(0.,10.))
include("plotting.jl")
plot_res(result, powergrid,H,disturbed_node)
