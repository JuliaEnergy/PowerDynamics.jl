using Test: @test
using LightGraphs: edges, Edge, SimpleGraph
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, ChangeInitialConditions, Inc, Dec, simulate, PowerGrid, State

graph = SimpleGraph(2)
Y = 0 + 5*im
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(nodes, lines)
u_Sl,u_Sw1,u_Sw2 = rand(ComplexF64, 3)
omega1 = 0.01
omega1_delta = 0.2
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), -omega1])

state2 = ChangeInitialConditions(2, :ω, Inc(omega1_delta))(state)
@test state2[2, :ω] == omega1 + omega1_delta

state3 = ChangeInitialConditions(2, :ω, Dec(omega1_delta))(state)
@test state3[2, :ω] == omega1 - omega1_delta

sol = simulate(ChangeInitialConditions(2, :ω, Inc(omega1_delta)),grid,state,timespan=(0,0.1))
@test sol !== nothing
@test sol.dqsol.retcode == :Success