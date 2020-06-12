using PowerDynamics
#using PowerDynamics: OperationPointError, find_operationpoint, rhs, RootRhs, SlackAlgebraic, SwingEqLVS, PQAlgebraic, StaticLine, PowerGrid, find_operationpoint, RootRhs, rhs, systemsize, find_operationpoint
using Test: @test, @testset, @test_throws, @test_logs, FallbackTestSetException

# define small test system
begin
    U1 = complex(1.0)
    P2 = -1.0
    Q2 = 0.0
    Y = 20f0im
    I_r = complex(0.5)
    H = 1.0
    P = 1.0
    D = 10.
    Ω = 2π*50
    Γ = 100.
    V = 1.

    nodes = [
        SlackAlgebraic(U = U1),
        PQAlgebraic(P = P2, Q = Q2),
        CSIMinimal(I_r = I_r),
        SwingEqLVS(H = H, P = P, D = D, Ω = Ω, Γ = Γ, V = V),
    ]
    lines = [
        StaticLine(from = 1, to = 2, Y = Y),
        StaticLine(from = 2, to = 3, Y = Y),
        StaticLine(from = 3, to = 4, Y = Y),
    ]

    grid = PowerGrid(nodes, lines)
end

@testset "type-specific initial condition guesses" begin
    @test_throws AssertionError PowerDynamics.guess.(grid.nodes, 5.0)

    expected_guess = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
    v = complex.(ones(length(nodes)))

    @test PowerDynamics.initial_guess(grid) == expected_guess
    @test PowerDynamics.initial_guess(grid, v) == expected_guess

    v[1] *= 0.9

    @test_throws AssertionError PowerDynamics.initial_guess(grid, v)

@testset "test for slack warning and SwingEq" begin
    no_slack = copy(nodes)
    no_slack[1] = PQAlgebraic(P=0., Q=0.)
    @test_logs (
        :warn,
        "There is no slack bus in the system to balance powers. Default voltage guess: u = 1 + 0j [pu].",
    ) PowerDynamics.initial_guess(PowerGrid(no_slack, lines));


    # missing slack should trigger warnings
    @test_logs (
            :warn,
            "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find an operation point.",
        )
        (
            :warn,
            "There is no slack bus in the system to balance powers. Default voltage guess: u = 1 + 0j [pu].",
        )
        find_operationpoint(PowerGrid(no_slack, lines));

    # test check for SwingEq nodes
    with_swing = copy(nodes)
    with_swing[end] = SwingEq(H = H, P = P, D = D, Ω = Ω)
    @test_throws OperationPointError find_operationpoint(PowerGrid(with_swing, lines))
end

@testset "check found operationpoint" begin

    # test default
    @test find_operationpoint(grid) == find_operationpoint(grid; method = :rootfind)

    # test case that should not converge
    should_fail = copy(nodes)
    should_fail[1] = PQAlgebraic(P=-10., Q=-10.)
    @test_throws OperationPointError find_operationpoint(PowerGrid(should_fail, lines); method = :nlsolve)

    op = find_operationpoint(grid; method = :nlsolve)

    @test op isa State

    # test if voltages are reasonable
    @test all(0.9 .< op[:, :v] .< 1.1)
    # frequency should be 0
    @test isapprox(op[4, :ω], 0., atol=1e-8)
    # Kirchhoff current law
    @test op[:, :i] |> sum == 0.
    # check set points
    @test isapprox(P2, op[2, :p], atol=1e-8)
    @test isapprox(Q2, op[2, :q], atol=1e-8)
    @test isapprox(U1, op[1, :u], atol=1e-8)
    @test isapprox(V, op[4, :v], atol=1e-8)

    op_rf = find_operationpoint(grid; method = :rootfind)

    # check that rootfind in the SteadyStateProblem is consistent with manual call to NLsolve
    @test op_rf.vec == op.vec

    # the following should not throw an OperationPointError since the SteadyStateProblem
    # still finds a solution where nlsolve fails in the current setup
    @test find_operationpoint(PowerGrid(should_fail, lines); method = :rootfind) isa State

    op_st = find_operationpoint(grid, PowerDynamics.initial_guess(grid); method = :steadystate)

    @test op_st.vec !== op.vec

    # frequency should be 0
    @test isapprox(op_st[4, :ω], 0., atol=1e-8)
    # Kirchhoff current law
    @test op_st[:, :i] |> sum == 0.

    # the voltage solution, however, is out of reasonable bounds
    @test any(0.9 .> op_st[:, :v]) | any(op_st[:, :v] .> 1.1)
end
