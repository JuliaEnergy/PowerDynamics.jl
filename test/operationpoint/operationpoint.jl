using PowerDynamics
#using PowerDynamics: OperationPointError, find_operationpoint, rhs, RootRhs, SlackAlgebraic, SwingEqLVS, PQAlgebraic, StaticLine, PowerGrid, find_operationpoint, RootRhs, rhs, systemsize, find_operationpoint
using LightGraphs: SimpleGraph, add_edge!
using Test: @test, @testset

U1 = complex(1., 0.)
P2 = -1.
Q2= 0.
Y = 20f0im
nodes = [SlackAlgebraic(U=U1), PQAlgebraic(P=P2,Q=Q2)]
lines = [StaticLine(from=1, to=2, Y=Y)]
graph = SimpleGraph(2)
add_edge!(graph, 1, 2)
grid = PowerGrid(graph, nodes, lines)
rpg = rhs(grid)
PowerDynamics.guess.(grid.nodes, 5.)
PowerDynamics.initial_guess(grid)
systemsize(grid)
op = find_operationpoint(grid)

op = find_operationpoint(grid; method = :nlsolve)

using OrdinaryDiffEq: ODEProblem, Rodas5
using SteadyStateDiffEq

op_prob = ODEProblem(rpg, PowerDynamics.initial_guess(grid), Inf)
sol = solve(
    SteadyStateProblem(op_prob),
    DynamicSS(Rodas5()),
)
sol.retcode == :Success

op2 = find_operationpoint(grid; method = :steadystate, abstol = 1e-6, reltol = 1e-6, tspan = Inf)

all(op[:, :v] .≈ op2[:, :v])
isapprox.(op[:, :φ], op2[:, :φ], atol=1e-8) |> all

op3 = find_operationpoint(grid; method = :rootfind)


@testset "algebraic only" begin
    U1 = 1+5im
    P2 = -2
    Q2= 3
    Y = 2+1.5im
    nodes = [SlackAlgebraic(U=U1), PQAlgebraic(P=P2,Q=Q2)]
    lines = [StaticLine(from=1, to=2, Y=Y)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2)
    grid = PowerGrid(graph, nodes, lines)
    op = find_operationpoint(grid)
    root = RootRhs(rhs(grid))
    @test all(root(convert(AbstractVector{Float64}, op)) .- zeros(systemsize(grid)) .< 1e-8)
    @test P2 ≈ op[2, :p]
    @test Q2 ≈ op[2, :q]
end

@testset "no convergence" begin
    P2 = -2
    Q2= 3im
    Y = 2+1.5im
    nodes = [PQAlgebraic(P=P2,Q=Q2), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=2)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2);
    lines = [StaticLine(from=1, to=2, Y=Y)]
    grid = PowerGrid(graph, nodes, lines)
    @test_throws OperationPointError find_operationpoint(grid)
end

@testset "algebraic and dynamic node" begin
    U1 = 1+5im
    P2 = -1
    V2 = 2
    nodes = [SlackAlgebraic(U=U1), SwingEqLVS(H=1, P=P2, D=1, Ω=50, Γ=20, V=V2)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2);
    lines = [StaticLine(from=1, to=2, Y=im)]
    grid = PowerGrid(graph, nodes, lines)
    op = find_operationpoint(grid)
    root = RootRhs(rhs(grid))
    @test all(root(convert(AbstractVector{Float64}, op)) .- zeros(systemsize(grid)) .< 1e-8)
    @test V2 ≈ op[2, :v]
    @test P2 ≈ op[2, :p]
end

@testset "error that SwingEqLVS should be used instead of SwingEq" begin
    nodes = [SwingEq(H=2, P =2, D=1, Ω=50), SlackAlgebraic(U=1)]
    graph = SimpleGraph(2)
    grid = PowerGrid(graph, nodes, [])

    @test_throws OperationPointError find_operationpoint(grid)
end
