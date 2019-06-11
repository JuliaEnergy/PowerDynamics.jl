using CSV
using DataFrames
using PowerDynBase
using SparseArrays
using NetworkDynamics
using DifferentialEquations

include("plotting.jl")

powergrid = read_network_from_csv("IEEE14_busses.csv", "IEEE14_lines.csv")

operationpoint = find_operationpoint(powergrid)

result = simulate(Perturbation(1, :ω, Inc(0.2)), powergrid, operationpoint, timespan = (0.0,0.3))
plot_res(result, powergrid)

result = simulate(LineFault(from=1,to=5), powergrid, operationpoint, timespan = (0.0,1.0))
plot_res(result, powergrid)
