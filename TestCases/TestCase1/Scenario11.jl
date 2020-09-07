using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise

perturbed_node=2

#grid frequency
ω = 2π*50

# per unit transformation base values
V_base_kV = 0.4
S_base_kW = 1
Y_base = S_base_kW*1000/(V_base_kV*1000)^2

# line paramterization
R_12 = 0.25
L_12 = 0.98*1e-3
X_12 = ω*L_12
Z_12 = R_12 + 1im*X_12
Y_12 = (1/Z_12)/Y_base

# node paramterization
P_2 = -16.67/S_base_kW
Q_2 = 0/S_base_kW

node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[VoltageDependentLoad(P=P_2,Q=Q_2,U=0.90,A=1.0,B=0.0)])
    # A = 0, B = 0 -> constant power load
    # A = 1, B = 0 -> constant impedance load
    # A = 0, B = 1 -> constant current load
line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=Y_12)])

powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)
startpoint = State(powergrid, [1.0, 0.0, 0.95, -0.04])

timespan = (0., 22.)

pd = PowerPerturbation(
    node = perturbed_node,
    fault_power = -11.11/S_base_kW,
    tspan_fault = (9.,18.),
    var = :P)

result_pd = simulate(pd, startpoint, timespan)

include("../../plotting.jl")
plot_res(result_pd,powergrid,perturbed_node)
