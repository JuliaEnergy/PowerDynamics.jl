using Pkg
Pkg.instantiate()
using PowerDynamics

node_list=[]

append!(node_list,[SlackAlgebraic(U=1.)])
append!(node_list,[PVInverterWithFrequencyControl(I_n=1.,k_PLL=0.1,f=50,f_s=50.01,T_m=0.2,k_P=0.2)])

line_list=[]
append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.2)])

powergrid = PowerGrid(node_list,line_list)

operationpoint = find_operationpoint(powergrid)

#result = simulate(Perturbation(2, :ω, Inc(0.5)), powergrid, operationpoint, timespan = (0.0,1.))
