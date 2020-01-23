using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise

perturbed_node=2

# per unit transformation base values
V_base_kV = 0.4
S_base_kW = 1
Y_base = S_base_kW/V_base_kV^2

# line paramterization
Y_12 = 1/((0.25+1im*0.98*1e-6)/Y_base)


node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[PQAlgebraic(P=-16.67/S_base_kW,Q=0./S_base_kW)])
line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=Y_12)])

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
