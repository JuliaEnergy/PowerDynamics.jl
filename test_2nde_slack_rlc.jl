#using Pkg
#Pkg.activate(".")
using PowerDynamics

admittance = 1im/0.2 #+ 0.03

node_list = []
append!(node_list, [SlackAlgebraic(U=1)])
append!(node_list, [RLCLoad(R=0.5,L=0.5,C=0.5)])


line_list =[StaticLine(Y=admittance, from=1, to=2)]
powergrid = PowerGrid(node_list, line_list)

operationpoint = find_operationpoint(powergrid)

result = simulate(Perturbation(1, :ω, Inc(0.2)), powergrid, operationpoint, timespan = (0.0,0.3))
plot_res(result, powergrid)
