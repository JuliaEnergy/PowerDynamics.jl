using Pkg
Pkg.instantiate()
using PowerDynamics
#using Revise
using OrderedCollections:OrderedDict

##

node_list=[]
    append!(node_list, [SlackAlgebraic(U=1.)])
    append!(node_list, [SwingEqLVS(H=5., P=1., D=0.1, Ω=50,Γ=0.1,V=1.)])
    append!(node_list, [SwingEqLVS(H=5., P=1., D=0.1, Ω=50,Γ=0.1,V=1.)])
    append!(node_list, [PQAlgebraic(P=-1,Q=-1)])

node_dict=OrderedDict(
    "bus1"=>SlackAlgebraic(U=1.),
    "bus2"=>SwingEqLVS(H=5., P=1., D=0.1, Ω=50,Γ=0.1,V=1.),
    "bus3"=>SwingEqLVS(H=5., P=1., D=0.1, Ω=50,Γ=0.1,V=1.),
    "bus4"=> PQAlgebraic(P=-1,Q=-1))

line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=-1im/0.02)])
    append!(line_list,[StaticLine(from=2,to=3,Y=-1im/0.02)])
    append!(line_list,[StaticLine(from=3,to=4,Y=-1im/0.02)])

line_dict=OrderedDict(
    "line1" => StaticLine(from="bus1",to="bus2",Y=-1im/0.02),
    "line2" => StaticLine(from="bus2",to="bus3",Y=-1im/0.02),
    "line3" => StaticLine(from="bus3",to="bus4",Y=-1im/0.02))


powergrid = PowerGrid(node_list,line_list)
powergrid_dict = PowerGrid(node_dict,line_dict)

##

operationpoint = find_operationpoint(powergrid)
operationpoint_dict = find_operationpoint(powergrid_dict)

pp = PowerPerturbation(node="bus2", fault_power=0.5,tspan_fault=(0.5,1.))
nsc = NodeShortCircuit(node="bus2", Y = complex(10., 0.),tspan_fault=(0.5,1.))
gp  = NodeParameterChange(node="bus2", value = 0.5,tspan_fault=(0.5,1.),var=:V)
lf=LineFailure(line_name="line1",tspan_fault=(0.5,1.))

##

#result = simulate(lf,operationpoint,(0.,2.))
data, result = power_flow(powergrid)

N = length(node_list)

v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:N]

@assert isapprox.(operationpoint[:, :v], v, atol=1e-10) |> all

va = [result["solution"]["bus"][string(k)]["va"] for k in 1:N]

@assert isapprox.(operationpoint[:, :φ], va, atol=1e-10) |> all

p = [result["solution"]["gen"][string(k)]["pg"] for k in 1:N-1]

@assert isapprox.(operationpoint[1:N-1, :p], p, atol=1e-10) |> all

q = [result["solution"]["gen"][string(k)]["qg"] for k in 1:N-1]

@assert isapprox.(operationpoint[1:N-1, :q], q, atol=1e-10) |> all

#result = simulate(LineFailure_new(1), powergrid, operationpoint, (0.,2.))
#result = simulate(lf, powergrid_dict, operationpoint, (0.,2.))
#result = simulate(pp_old,powergrid,operationpoint,(0.,2.))

#include("plotting_test.jl")
#plot_res(result,powergrid,2)
