using Test: @test
using PowerDynamics:
  simulate,
  PowerGrid,
  systemsize,
  rhs,
  CSIMinimal,
  SlackAlgebraic,
  StaticLine,
  NodeShortCircuit,
  find_valid_initial_condition
using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

busses = [
  SlackAlgebraic(; U = complex(1.0, 0.0)),
  CSIMinimal(; I_r = complex(0.4988, -0.4988)),
]


lines = [
  StaticLine(; from = 1, to = 2, Y = complex(1.736, -208.326)),
  PiModelLine(;
    from = 1,
    to = 2,
    y = complex(1.736, -208.326),
    y_shunt_km = 0.0,
    y_shunt_mk = 0.0,
  ),
]



nsc = NodeShortCircuit(node = 1, R = 0.1, sc_timespan=(1.,2.))

pg_static = PowerGrid(busses, stat_lines)
pg = PowerGrid(busses, lines)

ic_guess = ones(systemsize(pg)) + 0. * randn(systemsize(pg))
ic_guess[[3,6,9]] .= 0.
ic_guess = find_valid_initial_condition(pg,ic_guess)


s = simulate(nsc, pg, ic_guess, (0.,30.));
