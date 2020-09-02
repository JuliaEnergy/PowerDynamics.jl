using PowerDynamics
##
P, Q = rand(2)
U = 1
A = 0.
B = 0.
##
slack = SlackAlgebraic(U = U)
line = StaticLine(;from = 1, to = 2, Y = 10 / (0.1152 + im * 0.0458))
load = VoltageDependentLoad(P = P, Q = Q, U = U, A = A, B = B)
##
pg = PowerGrid([slack, load], [line,])
op = find_operationpoint(pg)
##
timespan = (0., 10.)
pd = PowerPerturbation(
    node = 2,
    fault_power = 0.,
    tspan_fault = (5.,6.),
    var = :P)

result_pd = simulate(pd, op, timespan)