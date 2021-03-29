using Test: @test, @testset
using PowerDynamics

# read grid data from json file
pg_read = read_powergrid("./examples/ieee14bus/grid.json", Json)

# build grid using the example script
include("../../examples/ieee14bus/buildexample.jl")
pg_build = powergrid

@testset "IEEE 14-bus stucture tests" begin
    @test pg_build.nodes == pg_read.nodes
    @test pg_build.lines == pg_read.lines
    @test pg_build.graph == pg_read.graph
end

operationpoint = find_operationpoint(pg_build)
timespan= (0.0,5.)

fault1 = ChangeInitialConditions(node="bus1", var=:ω, f=Inc(0.2))
fault2 = LineFailure(line_name="branch2", tspan_fault=(1.,5.))
fault3 = PowerPerturbation(node="bus5", fault_power=0.0, tspan_fault=(1.,5.), var=:P)

@testset "IEEE 14-bus solution tests" begin
    solution1 = simulate(fault1, pg_build, operationpoint, timespan)
    @test solution1.dqsol.retcode == :Success

    solution2 = simulate(fault2, pg_build, operationpoint, timespan)
    @test solution2.dqsol.retcode == :Success

    solution3 = simulate(fault3, pg_build, operationpoint, timespan)
    @test solution3.dqsol.retcode == :Success
end