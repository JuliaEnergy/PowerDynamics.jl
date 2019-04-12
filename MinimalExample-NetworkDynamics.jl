# note: please run before this script the following line
#]add git@github.com:FHell/NetworkDynamics.jl.git
#PkgDev.generate("PowerDynBase", "network_dynamics")
using PowerDynBase
#include("src/DynamicNodeMacro.jl")
#include("src/NodeDynamics/SwingEquation.jl")
using NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots

pyplot()

g = barabasi_albert(10,5)

begin
    #=
    Conventon is the following:
    v[1] + 1.j*v[2] is the complex voltage at a vertex.
    e[1] + 1.j*e[2] is the complex current at an edge.
    =#
end
struct complex_admittance_edge!
    admittance
end

function (cae::complex_admittance_edge!)(e,v_s,v_d,p,t)
    source_voltage = v_s[1] + v_s[2]*im
    destination_voltage = v_d[1] + v_d[2]*im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = cae.admittance * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end

begin
    # Making sure the function works as intended
    v1 = [1.,2.]
    v2 = [3.,4.]
    e = zeros(10)
    e1 = view(l, 5:6)
    cae = complex_admittance_edge!(3. + 5*im)
    cae(e1,v1,v2,nothing,0.)
    println(e1)
end

struct PQVertex
    P_complex
end
function (pq::PQVertex)(dv, v, e_s, e_d, p, t)
    current = total_current(e_s, e_d)
    voltage = v[1] + v[2] * im
    residual = pq.P_complex - voltage * conj(current)
    dv[1] = real(residual)
    dv[2] = imag(residual)
    nothing
end

struct complex_admittance_edge!
    admittance
end

function (cae::complex_admittance_edge!)(e,v_s,v_d,p,t)
    source_voltage = v_s[1] + v_s[2]*im
    destination_voltage = v_d[1] + v_d[2]*im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = cae.admittance * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end

swing_par = SwingEq(H=0.1, P=1, D=0.1, Ω=50)
swing_dyn = construct_node_dynamics(swing_par)

# Example PQ node:
pq_1 = StaticVertex(f! = PQVertex(randn() + randn()*im),
                 dim = 2)

using GraphPlot
gplot(g)

pq_list = [ODEVertex(f! = PQVertex(randn() + randn()*im),
                     dim = 2,
                     massmatrix = 0.,
                     sym = [:v_r, :v_i])
           for i in 1:5]

vertex_list = [construct_node_dynamics(SwingEq(H=abs(0.1*randn()), P=1, D=abs(0.1*randn()), Ω=50))
              for i in 1:5]

append!(vertex_list, pq_list)


edge_list = [StaticEdge(f! = complex_admittance_edge!(0.0 - 5.0im),
                        dim = 2)
             for e in edges(g)]

power_network_rhs = network_dynamics(vertex_list, edge_list, g)

begin
    x0 = rand(25)
    test_prob = ODEProblem(power_network_rhs,x0,(0.,50.))
end

#TODO: solve UndefVarError here, total_current not defined
test_sol = solve(test_prob, Rosenbrock23(autodiff=false), force_dtmin=true)

struct root_rhs
    rhs
    mm
end
function (rr::root_rhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mm * dx .- dx
end

rr = root_rhs(power_network_rhs, power_network_rhs.mass_matrix)

using NLsolve

nl_res = nlsolve(rr, x0)
ic = nl_res.zero
test_prob = ODEProblem(power_network_rhs,ic,(0.,50.))
test_sol = solve(test_prob, Rosenbrock23(autodiff=false))

test_sol
plot(test_sol)

plot(test_sol, vars = [s for s in power_network_rhs.syms if occursin("ω", string(s))])

plot(test_sol, vars=[:ω_1])
