using PowerDynamics: OperationPointError, find_operationpoint, rhs, RootRhs, SlackAlgebraic, SwingEqLVS, PQAlgebraic, StaticLine, 
PowerGrid, find_operationpoint, RootRhs, rhs, systemsize
using Test: @test, @testset


@testset "algebraic only" begin
    U1 = 1+5im
    S2 = -2+3im
    Y = 2+1.5im
    nodes = [SlackAlgebraic(U=U1), PQAlgebraic(S=S2)]
    graph = SimpleGraph(2)
    add_edge!(graph, 1, 2);
    lines = [StaticLine(from=1, to=2, Y=Y)]
    grid = PowerGrid(graph, nodes, lines)
    op = find_operationpoint(grid)
    root = RootRhs(rhs(grid))
    @test all(root(convert(AbstractVector{Float64}, op)) .- zeros(systemsize(grid)) .< 1e-8)
    @test S2 ≈ op[2, :s]
end

@testset "no convergence" begin
    S2 = -2+3im
    Y = 2+1.5im
    nodes = [PQAlgebraic(S=S2), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=2)]
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

@testset "error that slack bus is missing" begin
    nodes = [PQAlgebraic(S=2), PQAlgebraic(S=-2)]
    graph = SimpleGraph(2)
    grid = PowerGrid(graph, nodes, [])

    @test_throws OperationPointError find_operationpoint(grid)
end

@testset "error that SwingEqLVS should be used instead of SwingEq" begin
    nodes = [SwingEq(H=2, P =2, D=1, Ω=50), SlackAlgebraic(U=1)]
    graph = SimpleGraph(2)
    grid = PowerGrid(graph, nodes, [])

    @test_throws OperationPointError find_operationpoint(grid)
end
