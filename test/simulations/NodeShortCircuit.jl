using Test: @test
using PowerDynamics: simulate, PowerGrid, systemsize, rhs, PiModelLine, PQAlgebraic, VSIMinimal, StaticLine, NodeShortCircuit
using OrdinaryDiffEq: ODEProblem, Rodas4
using NetworkDynamics: find_valid_ic
import DiffEqBase: solve

busses = [VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  PQAlgebraic(S=-1.0-im*0.5),
  PQAlgebraic(S=-1.0-im*0.5),
  PQAlgebraic(S=-1.0-im*0.5)]
  #RLC_Load(R,L,C),RLC_Load(R,L,C),RLC_Load(R,L,C)]


lines = [
  PiModelLine(;from=1, to=4, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=5, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=3, to=6, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=3, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=1, to=2, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.)]

stat_lines = [
  StaticLine(;from=1, to=4, Y = 1/(0.1152 + im*0.0458)),
  StaticLine(;from=2, to=5, Y = 1/(0.1152 + im*0.0458)),
  StaticLine(;from=3, to=6, Y = 1/(0.1152 + im*0.0458)),
  StaticLine(;from=2, to=3, Y = 1/(0.1152 + im*0.0458)),
  StaticLine(;from=1, to=2, Y = 1/(0.1152 + im*0.0458))]

nsc = NodeShortCircuit(R = 0.1, node = 1, sc_timespan=(1.,2.))

pg_static = PowerGrid(busses, stat_lines)
pg = PowerGrid(busses, lines)

@test nsc(pg_static) == nothing

ic_guess = ones(systemsize(pg)) + 0. * randn(systemsize(pg))
ic_guess[[3,6,9]] .= 0.
ic_guess = find_valid_ic(rhs(pg),ic_guess)
problem = ODEProblem(rhs(pg),ic_guess,(0.,200.))
sol = solve(problem, Rodas4(autodiff=false))
op_point = sol[end]

s1, s2, s3 = simulate(nsc, pg, op_point, (0.,30.));
