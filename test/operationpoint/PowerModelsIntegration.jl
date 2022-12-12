using PowerDynamics: SlackAlgebraic, SwingEqLVS, PQAlgebraic, StaticLine, PowerGrid, find_operationpoint, VSIMinimal, power_flow
using OrderedCollections: OrderedDict
using Test: @test, @testset, @test_throws, @test_logs, @test_nowarn

node_dict = OrderedDict(
    "bus1" => SlackAlgebraic(U = 1.),
    "bus2" => SwingEqLVS(H = 5., P = 1., D = 0.1, Ω = 50, Γ = 0.1, V = 1.),
    "bus3" => SwingEqLVS(H = 5., P = 1., D = 0.1, Ω = 50, Γ = 0.1, V = 1.),# VSIMinimal(τ_P=1., τ_Q=2., K_P=3., K_Q=4., V_r=1, P=1., Q=0.),
    "bus4" => PQAlgebraic(P = -1, Q = -1)) # note - for PowerModels it is different default direction, multiply by -1

line_dict = OrderedDict(
    "line1" => StaticLine(from = "bus1", to = "bus2", Y = -1im / 0.02),
    "line2" => StaticLine(from = "bus2", to = "bus3", Y = -1im / 0.02),
    "line3" => StaticLine(from = "bus3", to = "bus4", Y = -1im / 0.02))

node_list = collect(values(node_dict))

line_list = []
push!(line_list, StaticLine(from = 1, to = 2, Y = -1im / 0.02))
push!(line_list, StaticLine(from = 2, to = 3, Y = -1im / 0.02))
push!(line_list, StaticLine(from = 3, to = 4, Y = -1im / 0.02))

powergrid_dict = PowerGrid(node_dict, line_dict)
powergrid = PowerGrid(node_list, line_list)

##

operationpoint = find_operationpoint(powergrid_dict; sol_method=:rootfind)
operationpoint_pf = find_operationpoint(powergrid_dict; sol_method=:rootfind,solve_powerflow=true)


data, result = power_flow(powergrid)

# The problem is in this line.
data_2, result_2 = power_flow(powergrid_dict)

@test all(keys(data_2) .== keys(data))

# compare bus data
for (k, v) in data["bus"]
    @test all(keys(v) .== keys(data_2["bus"][k]))
    @test all(values(v) .== values(data_2["bus"][k]))
end

for (k, v) in data["gen"]
    @test all(keys(v) .== keys(data_2["gen"][k]))
    @test all(values(v) .== values(data_2["gen"][k]))
end

for (k, v) in data["load"]
    @test all(keys(v) .== keys(data_2["load"][k]))
    @test all(values(v) .== values(data_2["load"][k]))
end

# compare branch data
for (k, v) in data["branch"]
    @test all(keys(v) .== keys(data_2["branch"][k]))
    @test all(values(v) .== values(data_2["branch"][k]))
end

N = length(node_list)

v =   [  result["solution"]["bus"][string(k)]["vm"] for k in 1:N]
v_2 = [result_2["solution"]["bus"][string(k)]["vm"] for k in 1:N]


@test isapprox.(operationpoint[:, :v], v, atol = 1e-6) |> all
@test isapprox.(operationpoint[:, :v], v_2, atol = 1e-6) |> all
@test isapprox.(operationpoint_pf[:, :v], v, atol = 1e-10) |> all
@test isapprox.(operationpoint_pf[:, :v], v_2, atol = 1e-10) |> all


va =   [  result["solution"]["bus"][string(k)]["va"] for k in 1:N]
va_2 = [result_2["solution"]["bus"][string(k)]["va"] for k in 1:N]

@test isapprox.(operationpoint[:, :φ], va, atol = 1e-7) |> all
@test isapprox.(operationpoint[:, :φ], va_2, atol = 1e-7) |> all
@test isapprox.(operationpoint_pf[:, :φ], va, atol = 1e-10) |> all
@test isapprox.(operationpoint_pf[:, :φ], va_2, atol = 1e-10) |> all

p = [result["solution"]["gen"][string(k)]["pg"] for k in 1:N-1]

@test isapprox.(operationpoint[:, :p][1:N - 1], p, atol = 1e-5) |> all

q = [result["solution"]["gen"][string(k)]["qg"] for k in 1:N - 1]

@test isapprox.(operationpoint[:, :q][1:N - 1], q, atol = 1e-4) |> all
