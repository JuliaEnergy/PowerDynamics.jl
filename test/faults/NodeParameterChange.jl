using Test: @test, @testset, @test_throws
using PowerDynamics: NodeParameterChange, simulate, SwingEqLVS, SlackAlgebraic, StaticLine, PowerGrid, State,systemsize, FieldUpdateError

##

@testset "test Power Perturbation with NodeParameterChange simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    npc = NodeParameterChange(node = 1, value = 0.9, tspan_fault = (0.1,1), var=:H)
    sol = simulate(npc, state, (0., 1.))
    @test sol !== nothing
    @test sol.dqsol.retcode == :Success
end

## 

@testset "test Power Perturbation with NodeParameterChange simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    npc = NodeParameterChange(node = 1, value = -0.9, tspan_fault = (0.1,1), var=:P)
    sol = simulate(npc, state, (0., 1.))
    @test sol !== nothing
    @test sol.dqsol.retcode == :Success
end



@testset "test PowerPerturbation simulation node without P should throw PowerPerturbationError" begin
    nodes = [SlackAlgebraic(U=-1*im), SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    npc = NodeParameterChange(node = 1, value = 0.9, tspan_fault = (0.1,1), var=:P)
    @test_throws FieldUpdateError simulate(npc, grid, state, (0., 1.))
end
