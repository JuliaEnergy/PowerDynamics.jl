using Test: @test,@testset, @test_throws
using PowerDynamics: PowerPerturbation, simulate, SwingEqLVS, SlackAlgebraic, StaticLine, PowerGrid, State,systemsize, FieldUpdateError
using SciMLBase: successful_retcode

@testset "test PowerPerturbation simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    pp = PowerPerturbation(node = 1, fault_power = -0.9, tspan_fault = (0.1,1))
    sol = simulate(pp, state, (0., 1.))
    @test sol !== nothing
@test successful_retcode(sol.dqsol)
end

@testset "test PowerPerturbation simulation P of type Int should also be accepted" begin
    invalidP = 0.9
    nodes = [SwingEqLVS(H=1., P=invalidP, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-invalidP, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    pp = PowerPerturbation(node = 1, fault_power = Int(1), tspan_fault = (0.1,1))
    sol = simulate(pp, state, (0., 1.))
    @test sol !== nothing
    @test successful_retcode(sol.dqsol)
end

@testset "test PowerPerturbation simulation node without P should throw NodePerturbationError" begin
    nodes = [SlackAlgebraic(U=-1*im), SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    pp = PowerPerturbation(node = 1, fault_power = 0.9, tspan_fault = (0.1, 1))
    @test_throws FieldUpdateError simulate(pp, state, (0., 1.))
end
