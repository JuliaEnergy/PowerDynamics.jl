using Pkg
Pkg.activate(".")
using PowerModels
using Ipopt
#using BenchmarkTools

# https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg200/
FILE = "test_powermodels/case3.m"
RAW_file = "test_powermodels/case3.raw"

# need a mapping from PD to a network_data dict and vice versa
network_data = parse_file(FILE)
network_raw_data = parse_file(RAW_file)

# use rectangular AC power flow in current/voltage formulation
model = instantiate_model(network_data, IVRPowerModel, build_pf_iv)
raw_model = instantiate_model(network_raw_data, IVRPowerModel, build_pf_iv)

# this prints all the equations and constraints
println(raw_model.model)

# run the actual power flow
pf = optimize_model!(raw_model, optimizer=Ipopt.Optimizer)

vr = Dict(name => data["vr"] for (name, data) in pf["solution"]["bus"])
vi = Dict(name => data["vi"] for (name, data) in pf["solution"]["bus"])

loss_ac =  Dict(name => data["pt"]+data["pf"] for (name, data) in pf["solution"]["branch"])

print_summary(pf["solution"])

# we can also access the admittance matrix
Y = calc_admittance_matrix(network_raw_data)
Y.matrix


# now build corrsponsing PowerDynamics model
using PowerDynamics
#using Revise

# see network_raw_data["bus"] and network_raw_data["gen"]
node_list=[]
    append!(node_list, [SlackAlgebraic(U=network_raw_data["gen"]["1"]["vg"])]) # bus 1 is of type 3 which is a slack bus?
    #append!(node_list, [SwingEqLVS(H=5., P=network_raw_data["gen"]["1"]["pg"], D=0.1, Ω=50,Γ=0.1,V=network_raw_data["gen"]["1"]["vg"])]) # is of type 2 which is a PV bus?
    #append!(node_list, [PQAlgebraic(P=network_raw_data["gen"]["3"]["pg"],Q=network_raw_data["gen"]["3"]["qg"])]) # is of type 2 which is a PV bus?
    append!(node_list, [SwingEqLVS(H=5., P=network_raw_data["gen"]["3"]["pg"], D=0.1, Ω=50,Γ=0.1,V=network_raw_data["gen"]["3"]["vg"])]) # is of type 2 which is a PV bus?
    #append!(node_list, [PQAlgebraic(P=network_raw_data["gen"]["5"]["pg"],Q=network_raw_data["gen"]["5"]["qg"])]) # is of type 2 which is a PV bus?
    append!(node_list, [SwingEqLVS(H=5., P=network_raw_data["gen"]["5"]["pg"], D=0.1, Ω=50,Γ=0.1,V=network_raw_data["gen"]["3"]["vg"])]) # is of type 2 which is a PV bus?

line_list=[]
    append!(line_list,[StaticLine(from=1,to=2,Y=Y.matrix[1,2])])
    append!(line_list,[StaticLine(from=2,to=3,Y=Y.matrix[2,3])])
    append!(line_list,[StaticLine(from=1,to=3,Y=Y.matrix[1,3])])


powergrid = PowerGrid(node_list,line_list)

operationpoint = find_operationpoint(powergrid)
#result = simulate(Perturbation(2, :ω, Inc(0.1)), powergrid, operationpoint, timespan = (0.0,10))

#include("plotting.jl")
#plot_res(result,powergrid,2)
