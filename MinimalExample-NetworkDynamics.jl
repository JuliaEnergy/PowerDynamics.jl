# note: please run before this script the following line
#]add git@github.com:FHell/NetworkDynamics.jl.git
#PkgDev.generate("PowerDynBase", "network_dynamics")
using PowerDynBase
using NetworkDynamics
using LightGraphs
using DifferentialEquations
using LinearAlgebra
# TODO: remove Plots package before merging into master
using Plots

pyplot()

g = barabasi_albert(10,5)


pq_list = [construct_node_dynamics(PQAlgebraic(S=1+im*1)) for i in 1:5]
vertex_list = [construct_node_dynamics(SwingEq(H=abs(0.1*randn()), P=1, D=abs(0.1*randn()), Ω=50))
              for i in 1:5]
append!(vertex_list, pq_list)

edge_list = [construct_edge(StaticLine(Y=0.0 - 5.0im)) for e in edges(g)]

power_network_rhs =
network_dynamics(vertex_list, edge_list, g)

begin
    x0 = rand(25)
    test_prob = ODEProblem(power_network_rhs,x0,(0.,50.))
end

test_sol = solve(test_prob, Rosenbrock23(autodiff=false), force_dtmin=true)

ic = find_valid_ic(power_network_rhs, x0)
test_prob = ODEProblem(power_network_rhs,ic,(0.,500.))
test_sol = solve(test_prob, Rosenbrock23(autodiff=false))

test_sol
plot(test_sol)

plot(test_sol, vars = [s for s in power_network_rhs.syms if occursin("ω", string(s))])

plot(test_sol, vars=[:ω_1])
