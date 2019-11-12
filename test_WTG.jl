using Pkg
Pkg.instantiate()
using PowerDynamics
using Revise


node_list=[]
    append!(node_list,[SlackAlgebraic(U=1.)])
    append!(node_list, [SwingEqLVS(H=5, P=-1., D=0.1, Ω=2π*50,Γ=0.1,V=1.)])
    #append!(node_list,[PQAlgebraic(P=-0.5,Q=0)])
    #append!(node_list,[FESS(u_dc_ref=1,D=1,C=1,J=0.018,ω_m_ref=0,k_PLL=1,f=50,L_g=0.001,L_m=0.01,P_n=2,Ψ_m=0.12,R_m=1,R_g=0.1,K_m1=5,K_m2=50,K_m3=0.1,K_m4=1,K_g1=5,K_g2=1000,K_g3=0.5,K_g4=1)])
    #append!(node_list,[PVInverterWithFrequencyControl(I_n=1.,k_PLL=20.,f=50,f_s=50.05,T_m=0.1,k_P=0.2,τ_ω=1.)])
    append!(node_list, [WindTurbineGenType4_RotorControl(T_L=0.0001,ω_0=2π*50,T_H=0.0055,K_ω=10,K_PPLL=50,K_IPLL=625,Q_ref=0,C=1,J=4,P=1,ω_rref=1,u_dcref=1,K_Q=0.1,K_v=400,K_pv=10,K_iv=10,K_ptrq=3*10,K_itrq=0.6*10)])
    #append!(node_list,[WindTurbineGenType4_RotorControl(T_L=1,K_dbr=0.0025,T_H=5.5,K_ω=10,K_PPLL=18.4,K_IPLL=170,Q_ref=0,C=1,J=4,P=1,ω_rref=1,u_dcref=1,K_Q=0.1,K_v=40,K_pv=10,K_iv=10,K_ptrq=3,K_itrq=0.6)])

    #append!(node_list, [WindTurbineGenType4_RotorControl(T_L=1,T_H=0.5,K_P=0.2,K_PLL=1,Q_ref=0,C=1,J=1,P=1,ω_rref=1,u_dcref=1,K_Q=1,K_v=1,K_pv=1,K_iv=100,K_pp=1,K_ip=10)])
    #append!(node_list, [WindTurbineGenType4(ΔP_max=0.3,k_P=0.2,K_PLL=10,Q_ref=0,C=1,J=1,P=1,ω_rref=1,K_Q=1,K_v=1,K_g1=1,K_g2=10,K_r1=1,K_r2=1)])
    line_list=[]
    append!(line_list,[PiModelLine(from=1, to=2, y=-1im/0.02, y_shunt_km=0.02*1im, y_shunt_mk=0.02*1im)])
    #append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.02)])
    append!(line_list,[PiModelLine(from=2, to=3, y=-1im/0.02, y_shunt_km=0.02*1im, y_shunt_mk=0.02*1im)])

powergrid = PowerGrid(node_list,line_list)


#using PowerDynamics: RootRhs
#using NLsolve: nlsolve, converged

#rr = RootRhs(rhs(powergrid))
#nl_res = nlsolve(rr, ic_guess)

operationpoint = find_operationpoint(powergrid)


disturbed_node=2;
final_time=2
simulation_time = (0.0,final_time);
pd = PowerPerturbation(node_number=disturbed_node,fraction=0.9, tspan_fault=(0.5,2));

#include("simulate_pd.jl")
include("plotting2.jl")
#sol_pd,result_pd = simulate_pd(pd,powergrid,operationpoint,simulation_time)
#plot_res(result_pd,powergrid,disturbed_node)
#result = simulate(pd,powergrid, operationpoint, simulation_time)
using OrdinaryDiffEq: ODEProblem, Rodas4, init, solve, step!, reinit!, savevalues!, u_modified!
p = Perturbation(disturbed_node, :ω, Dec(0.5));
x1 = p(operationpoint)
problem = ODEProblem{true}(rhs(powergrid),x1.vec,simulation_time)
sol = solve(problem,Rodas4(autodiff=false),saveat=1e-4)#,callback=cb)
result=PowerGridSolution(sol, powergrid)
plot_res(result,powergrid,disturbed_node,simulation_time)

#import DiffEqBase: solve
#Susing DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra
#saved_values = SavedValues(Float64, Vector{Float64})
#cb = SavingCallback((u,t,integrator)->(get_du(integrator), saved_values))


#p1 = plot(result, [2], :θ_PLL)
#    p2 = plot(result, [1], :ω)
#    p3 = plot(result2(range(0.,stop=3.,length=10000),Val{1},idxs=[6]))
    #p4 = plot(test[:,2])
#    p5 = plot(result2(range(0.,stop=3.,length=10000),Val{1},idxs=[3]))
#    #p6 = plot(test[:,3])
#    plot(
#    p1,p2,p3,p5;
#    layout=(4,1),
#    size = (500, 500))



#using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, ODEInterfaceDiffEq
#prob = ODEProblem{true}(rhs(powergrid), x1.vec, (0.0,1.0))
#saved_values = SavedValues(Float64, Array{Float64})
#cb = SavingCallback((u,t,integrator)-> integrator(t,Val{1}), saved_values,saveat=1e-6:1e-6:1.)
#sol = solve(prob, Ros4LStab(autodiff=false),saveat=1e-7,callback=cb)
#test=hcat(saved_values.saveval...)'
#plot(test[:,6])
#include("helpers.jl")
#r,n=determine_rocof_nadir(powergrid,sol,final_time)

#test=sol(range(0.,stop=final_time,length=10000),Val{1},idxs=variable_index(powergrid.nodes, 2, :θ_PLL)).u
#plot(result,[2],:p)
#plot(test)
