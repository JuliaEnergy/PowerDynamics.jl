using Pkg
Pkg.instantiate()
using PowerDynamics

node_list=[]

append!(node_list,[SlackAlgebraic(U=1.)])
append!(node_list, [SwingEqLVS(H=1, P=-1, D=0.01, Ω=50,Γ=0.1,V=1.)])
append!(node_list,[PVInverterWithFrequencyControl(I_n=1.,k_PLL=0.1,f=50,f_s=50,T_m=0.2,k_P=0.2)])

line_list=[]
append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.2)])
append!(line_list,[StaticLine(from=2,to=3,Y=-1im/0.2)])


powergrid = PowerGrid(node_list,line_list)

#vertices = map(construct_vertex, powergrid.nodes)
#println(v.mass_matrix for v in vertices)

operationpoint = find_operationpoint(powergrid)

#esult = simulate(Perturbation(2, :ω, Inc(0.5)), powergrid, operationpoint, timespan = (0.0,1.))

import DiffEqBase: solve
result = solve(powergrid, operationpoint,(0.0,1.))

include("plotting.jl")
plot_res(result,powergrid,2)

ω_indices = findall(n -> isa(n, PVInverterWithFrequencyControl), powergrid.nodes)
pl_ω = plot(result, ω_indices, :P)
