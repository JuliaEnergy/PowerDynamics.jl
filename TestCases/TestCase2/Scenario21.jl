using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise

perturbed_node=2

node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[PQAlgebraic(P=-0.5,Q=0.)])
    append!(node_list,[GridFormingTecnalia(τ_U, τ_I, τ_P, τ_Q, n_P, n_Q, k_P, k_Q, P, Q, V_r, R_f, X_)])
line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.02)])

powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)

timespan = (0., 20.)
pd = PowerPerturbation(
    fraction = 0.7,
    node_number = perturbed_node,
    tspan_fault = (1.,10.))

result_pd = simulate(pd,
    powergrid, operationpoint, timespan)

include("../../plotting.jl")
plot_res(result_pd,powergrid,perturbed_node)
