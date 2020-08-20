using Test: @test,@testset, @test_throws
using PowerDynamics: NodeParameterChange, simulate, SwingEqLVS, SlackAlgebraic, StaticLine, PowerGrid, State,systemsize, PowerPerturbationError

@testset "test Power Perturbation with GenericPertubation simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    sol = simulate(NodeParameterChange(
        node = 1,
        var_new = 0.9,
        tspan_fault = (0.1,1)),
        grid, state, (0., 1.))
    @test sol !== nothing
    @test sol.dqsol.retcode == :Success
end

@testset "test Power Perturbation with GenericPertubation simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    sol = simulate(NodeParameterChange(
        node = 1,
        var_new = 0.9,
        tspan_fault = (0.1,1),
        var=:P),
        grid, state, (0., 1.))
    @test sol !== nothing
    @test sol.dqsol.retcode == :Success
end



@testset "test PowerPerturbation simulation node without P should throw PowerPerturbationError" begin
    nodes = [SlackAlgebraic(U=-1*im), SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    @test_throws NodePerturbationError simulate(NodeParameterChange(
        node = 1,
        power_new = 0.9,
        tspan_fault = (0.1,1),
        var=:P),
        grid, state, (0., 1.))
end
