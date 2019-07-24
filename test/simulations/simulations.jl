using Test: @test
using LightGraphs: SimpleGraph, add_edge!, edges, Edge
using PowerDynamics: SlackAlgebraic, SwingEqLVS, StaticLine, Perturbation, Inc, Dec, LineFault

Y = 0 + 5*im
graph = SimpleGraph(3)
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]

add_edge!(graph, 1, 2)
add_edge!(graph, 1, 3)
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(graph, nodes, lines)
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

@testset "test SinglePhaseShortCircuitToGround simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2);
    lines = [StaticLine(Y=-im)]
    grid = PowerGrid(graph, nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    sol = simulate(SinglePhaseShortCircuitToGround(
        from=1,
        to=2,
        line_fraction=0.99,
        short_circuit_admittance=1e6,
        t_fault=.1,
        t_clearing=.5)
    @test sol != nothing
end

@testset "test PowerDrop simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2);
    lines = [StaticLine(Y=-im)]
    grid = PowerGrid(graph, nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    sol = simulate(PowerDrop(
        fraction = 0.9,
        node_number = 1,
        t_fault = 2.,
        t_clearing=3.),
        grid, state, (0., 5.))
    @test sol != nothing
end
