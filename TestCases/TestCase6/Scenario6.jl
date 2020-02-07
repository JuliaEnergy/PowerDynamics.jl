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
# line paramterization
Y_12 = (1/(0.25+1im*0.98*ω*1e-6))/Y_base
Y_13 = (1/(0.0474069+1im*0.0645069))/Y_base
Y_13_shunt = (1/((536.8837037+125.25700254999994*1im)))/Y_base

# node powers
P_2 = -16.67/S_base_kW
Q_2 = 0/S_base_kW
P_3=20/S_base_kW
Q_3=0/S_base_kW

# paramterization of grid-follwing inverter
ω_r = 0
ω_ini=0
w_cU=1.9998000 #rad: Voltage low pass filter cutoff frequency
τ_U = 1/w_cU
V_r =1 # pu
P_r=20 # pu
Q_r=0 # pu

K_pω=2π*0.001
K_iω=2π*0.02
K_ω=1/(2π*0.1/(40000))
K_v= 1/(0.916*398*(1.1-0.9)/(40000- (-40000) ))


node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[PQAlgebraic(P=P_2,Q=Q_2)])
    append!(node_list,[GridFollowingTecnalia(τ_u=τ_U,ω_ini=ω_ini,K_pω=K_pω,K_iω=K_iω,K_ω=K_ω,K_v=K_v,ω_r=ω_r,V_r=V_r,P=P_r,Q_r=Q_r)])
    append!(node_list,[PQAlgebraic(P=-0.1,Q=-0.1)])
    #append!(node_list,[Connector()])

line_list=[]
    append!(line_list,[PiModelLine(from=3,to=4,y=Y_13, y_shunt_km=Y_13_shunt, y_shunt_mk=Y_13_shunt)])
    append!(line_list,[StaticLine(from=2,to=4,Y=Y_12)])
    append!(line_list,[StaticLine(from=1,to=4,Y=Y_13)])
    #append!(line_list,[ConnectorLine(from=1,to=4)])


powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)

timespan = (0., 20.)
pd = PowerPerturbation(
    fraction = 11.11/16.67,
    node_number = perturbed_node,
    tspan_fault = (1.,10.))

result_pd = simulate(pd,
    powergrid, operationpoint, timespan)

include("../../plotting.jl")
plot_res(result_pd,powergrid,perturbed_node)
