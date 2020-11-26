using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise

perturbed_node="VS Inverter"

#grid frequency
ω = 2π*50

# per unit transformation base values
V_base_kV = 0.4
S_base_kW = 1
Y_base = S_base_kW*1000/(V_base_kV*1000)^2

# line paramterization
R_14 = 0.25
L_14 = 0.98*1e-3
X_14 = ω*L_14
Z_14 = R_14 + 1im*X_14
Y_14 = (1/Z_14)/Y_base

# transformer parametrization
Y_34 = (1/(0.0474069+1im*0.0645069))/Y_base
Y_34_shunt = (1/((536.8837037+125.25700254999994*1im)))/Y_base

# node powers
P_2 = -16.67/S_base_kW
Q_2 = 0.0/S_base_kW
P_3 = 20/S_base_kW
Q_3 = 0/S_base_kW

# paramterization of grid-forming inverter
ω_r = 0#-0.25
w_cU=198.019802 #rad: Voltage low pass filter cutoff frequency
τ_U = 1/w_cU
w_cI=198.019802 #rad: Current low pass filter cutoff frequency
τ_I=1/w_cI
w_cV = 0.999000999
τ_V = 1/w_cV
w_cω = 0.999000999
τ_ω = 1/w_cV
w_cP=0.999000999 #rad: Active power low pass filter cutoff
τ_P=1/w_cP
w_cQ=0.999000999 #rad: Reactive power low pass filter cutoff
τ_Q =1/w_cQ
n_P=10 # Inverse of the active power low pass filter
n_Q=10 #Inverse of the reactive power low pass filter
K_P= 2π*(50.5-49.5)/(40/S_base_kW) # P-f droop constant Hz/kW is transformed into Hz/kW
transformer_ratio = 0.916
K_Q=1/transformer_ratio*(398/(1000*V_base_kV)*(1.1-0.9)/(40/S_base_kW- (-40/S_base_kW))) # Q-U droop constant, kV/kVar is transformed into pu/pu
R_f=0.6 #in Ω #Vitual resistance
X_f=0.8 # in Ω #Vitual reactance
V_r =1

node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[VoltageDependentLoad(P=P_2,Q=Q_2,U=1.,A=1.0,B=0.0)])
    append!(node_list,[GridFormingTecnalia(ω_r=0,τ_U=τ_U, τ_I=τ_I, τ_P=τ_P, τ_Q=τ_Q, n_P=n_P, n_Q=n_Q, K_P=K_P, K_Q=K_Q, P=P_3, Q=Q_3, V_r=V_r, R_f=R_f, X_f=X_f)])
    #append!(node_list,[VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P_3,Q=Q_3)])
    #append!(node_list,[Connector()])
line_list=[]
    #append!(line_list,[ConnectorLine(from=2,to=4)])
    append!(line_list,[PiModelLine(from=3,to=2,y=Y_34,y_shunt_mk=Y_34_shunt/2,y_shunt_km=Y_34_shunt/2)])
    append!(line_list,[StaticLine(from=1,to=2,Y=Y_14)])

node_dict=OrderedDict(
    "Slack"=> SlackAlgebraic(U=1.),
    "Load"=>VoltageDependentLoad(P=P_2,Q=Q_2,U=1.,A=1.0,B=0.0),
    "VS Inverter"=>GridFormingTecnalia(ω_r=0,τ_U=τ_U, τ_I=τ_I, τ_P=τ_P, τ_Q=τ_Q, n_P=n_P, n_Q=n_Q, K_P=K_P, K_Q=K_Q, P=P_3, Q=Q_3, V_r=V_r, R_f=R_f, X_f=X_f))
line_dict=OrderedDict(
    "Line1"=>StaticLine(from="Slack",to="Load",Y=Y_12),
    "Line2"=>PiModelLine(from="Load",to="VS Inverter",y=Y_23,y_shunt_mk=Y_23_shunt/2,y_shunt_km=Y_23_shunt/2))
    #"Line3"=>PiModelLine(from="Load",to="CS Inverter",y=Y_24,y_shunt_mk=Y_24_shunt/2,y_shunt_km=Y_24_shunt/2))


powergrid = PowerGrid(node_dict,line_dict)
operationpoint = find_operationpoint(powergrid)

timespan = (0., 40.)
P3_new = P_3-0.25/K_P*2*pi

pd = PowerPerturbation(
    node = perturbed_node,
    fault_power = P3_new/S_base_kW,
    tspan_fault = (12.,25.),
    var = :P)

result_pd = simulate(pd,
    powergrid, operationpoint, timespan)

include("../../plotting.jl")
create_plot(result_pd)
#plot_res(result_pd,powergrid,perturbed_node)
