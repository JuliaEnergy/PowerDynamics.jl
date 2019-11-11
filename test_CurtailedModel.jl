using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise


node_list=[]
    #append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list, [SwingEqLVS(H=0.5, P=-0.5, D=0.1, Ω=50,Γ=0.1,V=1.)])
    #append!(node_list,[PQAlgebraic(P=-1,Q=0)])
    #append!(node_list,[FESS(u_dc_ref=1,D=1,C=1,J=0.018,ω_m_ref=0,k_PLL=1,f=50,L_g=0.001,L_m=0.01,P_n=2,Ψ_m=0.12,R_m=1,R_g=0.1,K_m1=5,K_m2=50,K_m3=0.1,K_m4=1,K_g1=5,K_g2=1000,K_g3=0.5,K_g4=1)])
    #append!(node_list,[PVInverterWithFrequencyControl(I_n=1.,k_PLL=20.,f=50,f_s=50.05,T_m=0.1,k_P=0.2,τ_ω=1.)])
    append!(node_list, [CurtailedPowerPlantWithInertia(PICFunction=true,P=0.5,ω_0=2π*50,T_AI=10,K_PPLL=50,K_IPLL=625,T_d=0.3,T_f=0.4,K_PV=0.28,K_IV=186)])
    #append!(node_list, [CurtailedPowerPlantWithInertia(P=1,ω_0=2π*50,T_AI=10,K_PPLL=600,K_IPLL=50,T_d=0.38,T_f=0.5,K_PV=10,K_IV=10)])
    #append!(node_list, [CurtailedPowerPlantWithInertia(P=1,ω_0=2π*50,T_AI=10,K_PPLL=100,K_IPLL=100,T_d=0.0005,T_f=0.001,K_PV=1,K_IV=1)])
    #append!(node_list, [WindTurbineGenType4_RotorControl(T_L=1,T_H=0.5,K_P=0.2,K_PLL=1,Q_ref=0,C=1,J=1,P=1,ω_rref=1,u_dcref=1,K_Q=1,K_v=1,K_pv=1,K_iv=100,K_pp=1,K_ip=10)])
    #append!(node_list, [WindTurbineGenType4(ΔP_max=0.3,k_P=0.2,K_PLL=10,Q_ref=0,C=1,J=1,P=1,ω_rref=1,K_Q=1,K_v=1,K_g1=1,K_g2=10,K_r1=1,K_r2=1)])
    line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.02)])
    #append!(line_list,[StaticLine(from=2,to=3,Y=-1im/0.02)])

powergrid = PowerGrid(node_list,line_list)


operationpoint = find_operationpoint(powergrid)


disturbed_node=1;
final_time=10
simulation_time = (0.0,final_time);
pd = PowerPerturbation(node_number=disturbed_node,fraction=0.9, tspan_fault=(0.5,5));

using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve
#using Sundials: CVODE_BDF
p = Perturbation(disturbed_node, :ω, Dec(0.5));
x1 = p(operationpoint)
problem = ODEProblem{true}(rhs(powergrid),x1.vec,simulation_time)
sol = solve(problem,Rodas4(autodiff=false),saveat=1e-5)#,callback=cb)
result=PowerGridSolution(sol, powergrid);
include("plotting2.jl")
plot_res(result,powergrid,disturbed_node,simulation_time)

#include("helpers.jl")
#r,n=determine_rocof_nadir(powergrid,sol,final_time)

#test=sol(range(0.,stop=final_time,length=10000),Val{1},idxs=variable_index(powergrid.nodes, 2, :θ_PLL)).u
#plot(result,[2],:p)
#plot(test)
