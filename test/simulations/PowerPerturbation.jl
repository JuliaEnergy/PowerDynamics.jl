using Test: @test,@testset
using PowerDynamics: PowerPerturbation, simulate, SwingEqLVS, StaticLine, PowerGrid, State

@testset "test PowerPerturbation simulation" begin
    nodes = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    sol = simulate(PowerPerturbation(
        fraction = 0.9,
        node_number = 1,
        tspan_fault = (2.,3.)),
        grid, state, (0., 5.))
    @test sol != nothing
    @test sol.dqsol.retcode == :Success
end
