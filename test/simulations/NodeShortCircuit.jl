using Test: @test, @testset
using PowerDynamics:
  simulate,
  PowerGrid,
  CSIMinimal,
  SlackAlgebraic,
  PowerGrid,
  StaticLine,
  PiModelLine,
  Transformer,
  NodeShortCircuit,
  find_operationpoint
using OrdinaryDiffEq: ODEProblem, Rodas4
import DiffEqBase: solve

nodes = [
  SlackAlgebraic(; U = complex(1.0, 0.0)),
  CSIMinimal(; I_r = complex(0.4988, -0.4988)),
]

y = complex(1.736, -208.326)

SL = StaticLine(; from = 1, to = 2, Y = y)
PL = PiModelLine(;
    from = 1,
    to = 2,
    y = y,
    y_shunt_km = 0.0,
    y_shunt_mk = 0.0,
  )

@testset "test NodeShortCircuit construction" begin
    pg = PowerGrid(nodes, [SL,])

    nsc = NodeShortCircuit(;
        node_number = 2,
        Y = complex(160., 0.),
        tspan_fault = (0.5, 0.65),
    )

    # test if nsc returns the correct type
    @test nsc(pg) isa PowerGrid

    fault_pg = nsc(pg)

    # test if the shunt field is appropriately adjusted
    @test getfield(fault_pg.nodes[nsc.node_number], nsc.shunt_symbol) == complex(160., 0.)

    # test with PiModelLine or Transformer
    @test nsc(PowerGrid(nodes, [PL,])) isa PowerGrid

end

@testset "test NodeShortCircuit simulation" begin
    pg = PowerGrid(nodes, [SL,])
    op = find_operationpoint(pg)

    # make sure we obtain the right operating point
    @test op[2, :v] == 1.0024169204096718

    # a short circuit at the slack should have no effect
    nsc = NodeShortCircuit(;
        node_number = 1,
        Y = complex(160., 0.),
        tspan_fault = (0.5, 0.65),
    )

    sol = simulate(nsc, op, (0., 1.));

    @test sol != nothing
    @test sol.dqsol.retcode == :Success

    # no voltage drop
    @test all(sol(0.6, 1:2, :v) .≈ op[1:2, :v])

    # a short circuit at the non-slack node
    nsc = NodeShortCircuit(;
        node_number = 2,
        Y = complex(160., 0.),
        tspan_fault = (0.5, 0.65),
    )

    sol = simulate(nsc, op, (0., 1.))

    # observe voltage drop
    @test sol(0.6, 1, :v) == op[1, :v]
    @test sol(0.6, 2, :v) == 0.7918311818509268

    # test whether currents are determined correctly
    normal_current = op[2, :i]
    fault_current = pg.lines[1].Y * (sol(0.6, 2, :u) - sol(0.6, 1, :u))
    shunt_current = - nsc.Y * sol(0.6, 2, :u)
    @test abs(fault_current) > abs(normal_current)
    @test fault_current - shunt_current ≈ normal_current
end
