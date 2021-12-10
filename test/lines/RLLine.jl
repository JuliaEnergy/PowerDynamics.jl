using Test: @testset, @test, @test_throws, @test_broken
using PowerDynamics: RLLine, StaticLine, construct_edge, symbolsof, dimension
using NetworkDynamics: ODEEdge
using LinearAlgebra: I

include("LineTestBase.jl")



@testset "RLLine construct_edge" begin
        R = rand_real()
        L = rand_real()
        omega = rand_real()
        R_INVALID, L_INVALID, omega_INVALID = rand_negative(3)

        line = RLLine(from = 1, to = 2, R = R, L = L, ω0 = omega)
        edge = construct_edge(line)

        @test_throws AssertionError construct_edge(RLLine(
                from = 1,
                to = 2,
                R = R_INVALID,
                L = L,
                ω0 = omega,
        ))
        @test_throws AssertionError construct_edge(RLLine(
                from = 1,
                to = 2,
                R = R,
                L = L_INVALID,
                ω0 = omega,
        ))
        @test_throws AssertionError construct_edge(RLLine(
                from = 1,
                to = 2,
                R = R,
                L = L,
                ω0 = omega_INVALID,
        ))

        @test isa(edge, ODEEdge)
        @test symbolsof(line) == [:id, :iq, :id_r, :iq_r]
        @test dimension(line) == 4
        @test edge.mass_matrix == I

    # assure function call does not explode!
        smoketest_rhs(edge, int_x = [], int_dx = [])
end



@testset "RLLine should have the right fixed point" begin
        line = RLLine(from = 1, to = 2, R = 0.01, L = 0.1, ω0 = 100π)
        edge = construct_edge(line)

        # construct steady state
        #Z = [line.R -line.ω0*line.L; line.ω0*line.L line.R]
        Zinv = [line.R line.ω0 * line.L; -line.ω0 * line.L line.R] ./
               (line.R^2 + line.ω0^2 * line.L^2)
        v_s = [10.0; 5.0]
        v_d = [12.0; 2.0]

        e_star = -Zinv * (v_s .- v_d)

        e = [e_star; e_star]
        de = similar(e)
        edge.f(de, e, v_s, v_d, 0, 0)

        @test de[1] == 0
        @test de[2] == 0
        @test de[3] == 0
        @test de[4] == 0
end


@testset "The steady state of RLLine should coincide with StaticLine" begin
        line = RLLine(from = 1, to = 2, R = 0.01, L = 0.1, ω0 = 100π)
        rll_edge = construct_edge(line)

        sl = StaticLine(
                ;
                from = line.from,
                to = line.to,
                Y = inv(complex(line.R, line.ω0 * line.L)),
        )
        sl_edge = construct_edge(sl)

        v_s = [10.0; 5.0]
        v_d = [12.0; 2.0]

        arr = zeros(4)
        sl_edge.f(arr, v_s, v_d, 0, 0)

        # recapitulate RLLine steady state
        Zinv = [line.R line.ω0 * line.L; -line.ω0 * line.L line.R] ./
               (line.R^2 + line.ω0^2 * line.L^2)
        e_star = -Zinv * (v_s .- v_d)
        e = [e_star; -e_star]

        @test arr ≈ e
end

using PowerDynamics: PowerGrid,
                     dimension,
                     rhs,
                     State,
                     simulate,
                     SlackAlgebraic,
                     SwingEqLVS,
                     PowerPerturbation
using NetworkDynamics: find_valid_ic

# extend to account for line variables
import PowerDynamics: systemsize, dimension, AbstractLine

@testset "Do static and dynamics solutions agree for Slack & Swing" begin

        # "Simulate power drop"

        line = RLLine(from = 1, to = 2, R = 0.1, L = 0.01, ω0 = 100π)
        sl = StaticLine(
                from = line.from,
                to = line.to,
                Y = inv(complex(line.R, line.ω0 * line.L)),)

        busses = [
                SlackAlgebraic(U = complex(1.0, 0.0)),
                SwingEqLVS(H = 1.0, P = 1.0, D = 1.0, Ω = 100π, Γ = 10., V = 1.0),]
        pg = PowerGrid(busses, [line,])
        pg_sl = PowerGrid(busses, [sl,])

        v_s = [1.0; 0.0]
        v_d = [0.9; 0.5]
        Zinv = [line.R line.ω0 * line.L; -line.ω0 * line.L line.R] ./
               (line.R^2 + line.ω0^2 * line.L^2)
        e_star = -Zinv * (v_s .- v_d)

        ic_guess = [
                v_s[1],
                v_s[2],
                v_d[1],
                v_d[2],
                0.0,
                e_star[1],
                e_star[2],
                -e_star[1], # make sure initial conditions have correct sign
                -e_star[2],]
        ic = find_valid_ic(rhs(pg_sl), ic_guess[1:5])

        sol_sl = simulate(
                PowerPerturbation(
                        fault_power = 0.5,
                        node = 2,
                        tspan_fault = (0.1, 0.2),
                ),
                pg_sl,
                State(pg_sl, ic),
                (0.0, 1),)

        ic = find_valid_ic(rhs(pg), ic_guess)

        sol = simulate(
                PowerPerturbation(
                        fault_power = 0.5,
                        node = 2,
                        tspan_fault = (0.1, 0.2),
                ),
                pg,
                State(pg, ic),
                (0.0, 1),)
        for t = 0.1:0.1:1.
                @test all(sol.dqsol(t)[1:5] .- sol_sl.dqsol(t) .< 10^-3)
        end
end

# using Plots
# begin
#         p = plot(sol.dqsol, vars = 1:4)
#         plot!(sol.dqsol, vars = 6:7, c = :red, linestyle = :solid)
#         plot!(sol.dqsol, vars = 8:9, c = :blue, linestyle = :dash)
#         plot!(sol_sl.dqsol, vars =  1:4, c = "black", linestyle = :dash)
#         display(p)
# end

### Two swings

@testset "Do static and dynamics solutions agree for double Swing" begin
        # "Simulate power drop"

        line = RLLine(from = 1, to = 2, R = 0.1, L = 0.01, ω0 = 100π)
        sl = StaticLine(
                from = line.from,
                to = line.to,
                Y = inv(complex(line.R, line.ω0 * line.L)),)

        busses = [
                SwingEqLVS(H = 1.0, P = 1.0, D = 1.0, Ω = 100π, Γ = 10., V = 1.0),
                SwingEqLVS(H = 1.0, P = 1.0, D = 1.0, Ω = 100π, Γ = 10., V = 1.0),]
        pg = PowerGrid(busses, [line,])
        pg_sl = PowerGrid(busses, [sl,])

        v_s = [1.0; 0.0]
        v_d = [0.9; 0.5]
        Zinv = [line.R line.ω0 * line.L; -line.ω0 * line.L line.R] ./
               (line.R^2 + line.ω0^2 * line.L^2)
        e_star = -Zinv * (v_s .- v_d)

        ic_guess = [
                v_s[1],
                v_s[2],
                0.0,
                v_d[1],
                v_d[2],
                0.0,
                e_star[1],
                e_star[2],
                -e_star[1], # make sure initial conditions have correct sign
                -e_star[2],]

        ic = find_valid_ic(rhs(pg_sl), ic_guess[1:6])

        sol_sl = simulate(
                PowerPerturbation(
                        fault_power = 0.5,
                        node = 2,
                        tspan_fault = (0.1, 0.2),
                ),
                pg_sl,
                State(pg_sl, ic),
                (0.0, 1),)

        ic = find_valid_ic(rhs(pg), ic_guess)

        sol = simulate(
                PowerPerturbation(
                        fault_power = 0.5,
                        node = 2,
                        tspan_fault = (0.1, 0.2),
                ),
                pg,
                State(pg, ic),
                (0.0, 1),)

        for t = 0.1:0.1:1.
                @test all(sol.dqsol(t)[1:6] .- sol_sl.dqsol(t) .< 10^-3)
        end
end

# begin
#         p = plot(sol.dqsol, vars = [1,2,4,5])
#         plot!(sol.dqsol, vars = 7:8, c = :red, linestyle = :solid)
#         plot!(sol.dqsol, vars = 9:10, c = :blue, linestyle = :dash)
#         plot!(sol_sl.dqsol, vars =  [1,2,4,5], c = "black", linestyle = :dash)
#         display(p)
# end
