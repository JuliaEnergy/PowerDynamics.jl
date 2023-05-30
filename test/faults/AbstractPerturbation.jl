using Test: @test
using Graphs: edges, Edge, SimpleGraph
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, ChangeInitialConditions, Inc, Dec, simulate, PowerGrid, State, AbstractPerturbation

graph = SimpleGraph(2)
Y = 0 + 5 * im
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(nodes, lines)
u_Sl, u_Sw1, u_Sw2 = rand(ComplexF64, 3)
omega1 = 0.01
omega1_delta = 0.2
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), -omega1])

Base.@kwdef struct DummyPerturbation <: AbstractPerturbation
    tspan_fault
end

function (dp::DummyPerturbation)(powergrid)
    powergrid
end

sol = simulate(DummyPerturbation(tspan_fault=(0.0, 0.1)), grid, state, (0, 0.1); reltol=1e-6)
@test sol !== nothing
@test SciMLBase.successful_retcode(sol.dqsol.retcode)

sol = simulate(DummyPerturbation(tspan_fault=(0.0, 0.1)), state, (0, 0.1); reltol=1e-6)
@test sol !== nothing
@test SciMLBase.successful_retcode(sol.dqsol.retcode)

