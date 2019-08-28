using Test: @test
using LightGraphs: edges, Edge
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, Perturbation, Inc, Dec, LineFault, PowerPerturbation, simulate

Y = 0 + 5*im
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(nodes, lines)
@syms u_Sl u_Sw1 u_Sw2
omega1 = 0.01
omega1_delta = 0.2
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])

state2 = Perturbation(2, :ω, Inc(omega1_delta))(state)
@test state2[2, :ω] == omega1 + omega1_delta

state3 = Perturbation(2, :ω, Dec(omega1_delta))(state)
@test state3[2, :ω] == omega1 - omega1_delta

line_fault = LineFault(from=1, to=2)
faulty_grid = line_fault(grid)

@test length(faulty_grid.lines) == 1
@test collect(edges(faulty_grid.graph)) == [Edge(1,3)]
