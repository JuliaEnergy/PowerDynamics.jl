using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise

perturbed_node=3


#grid frequency
ω = 2π*50

# per unit transformation base values
V_base_kV = 0.4
S_base_kW = 1
Y_base = S_base_kW*1000/(V_base_kV*1000)^2

# line paramterization
Y_12 = (1/(0.25+1im*0.98*ω*1e-6))/Y_base
Y_13 = (1/(0.0474069+1im*0.0645069))/Y_base
# node paramterization
P_2 = -16.67/S_base_kW
Q_2 = 0/S_base_kW

# paramterization of grid-forming inverter

w_cU=198.019802 #Hz: Voltage low pass filter cutoff frequency
τ_U = 1/w_cU
w_cI=198.019802 #Hz: Current low pass filter cutoff frequency
τ_I=1/w_cI
w_cP=0.999000999 #Hz: Active power low pass filter cutoff
τ_P=1/w_cP
w_cQ=0.999000999 #Hz: Reactive power low pass filter cutoff
τ_Q =1/w_cQ
n_P=10 # Inverse of the active power low pass filter
n_Q=10 #Inverse of the reactive power low pass filter
K_P=2π*(50.5-49.5)/(40000/1000) # P-f droop constant
transformer_ratio = 1.08686
K_Q=0.916*398*(1.1-0.9)/(40000/1000- (-40000/1000)) # Q-U droop constant
R_f=0.6*ω #Vitual resistance
X_f=0.8*ω #Vitual reactance
V_r =1
P=20/S_base_kW
Q=0/S_base_kW

node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list,[PQAlgebraic(P=P_2,Q=Q_2)])
    append!(node_list,[GridFormingTecnalia_experimental(ω_r=0,τ_U=τ_U, τ_I=τ_I, τ_P=τ_P, τ_Q=τ_Q, n_P=n_P, n_Q=n_Q, K_P=K_P, K_Q=K_Q, P=P, Q=Q, V_r=V_r, R_f=R_f, X_f=X_f)])
    #append!(node_list,[VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)])
    append!(node_list,[Connector()])
line_list=[]
    append!(line_list,[ConnectorLine(from=1,to=4)])
    append!(line_list,[StaticLine(from=3,to=4,Y=Y_13)])
    append!(line_list,[StaticLine(from=2,to=4,Y=Y_12)])


powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)

timespan = (0., 15.)
pd = PowerPerturbation(
    node = perturbed_node,
    fault_power = 10/S_base_kW,
    tspan_fault = (1.,10.),
    var = :P)

result_pd = simulate(pd,
    powergrid, operationpoint, timespan)

include("../../plotting.jl")
plot_res(result_pd,powergrid,perturbed_node)
