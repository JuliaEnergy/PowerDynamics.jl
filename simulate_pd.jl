using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve!, step!, reinit!, savevalues!, u_modified!


function simulate_pd(pd,powergrid,operationpoint,timespan)
    x0 = operationpoint;
    problem = ODEProblem{true}(rhs(powergrid), x0.vec, timespan)
    integrator = init(problem, Rodas4(autodiff=false))

    step!(integrator, pd.tspan_fault[1], true)

    # update integrator with error
    integrator.f = rhs(pd(powergrid))
    u_modified!(integrator,true)

    step!(integrator, pd.tspan_fault[2]-pd.tspan_fault[1], true)

    # update integrator, clear error
    integrator.f = rhs(powergrid)
    u_modified!(integrator,true)

    step!(integrator, timespan[2]-pd.tspan_fault[2], true)


    solve!(integrator)
    result_pd2 = PowerGridSolution(integrator.sol, powergrid)
    return (integrator.sol,result_pd2)
end
