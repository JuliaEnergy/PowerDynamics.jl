using PowerDynamics: read_powergrid, write_powergrid, Json, PowerGrid, NormalForm
using OrderedCollections: OrderedDict
using Test: @test, @test_throws

expected_swing_eq = SwingEq(H = 1, P = 2, D = 3, Ω = 12)
expected_swing_eq_lvs = SwingEqLVS(H = 1, P = 2, D = 3, Ω = 12, Γ = 13, V = 14)
expected_fourth_order_eq = FourthOrderEq(
    H = 1,
    P = 2,
    D = 3,
    Ω = 12,
    E_f = 13,
    T_d_dash = 14,
    T_q_dash = 15,
    X_q_dash = 16,
    X_d_dash = 17,
    X_d = 18,
    X_q = 19,
)
expected_slack = SlackAlgebraic(U = 10)
expected_pq = PQAlgebraic(P = 0, Q = 10)
expected_pv = PVAlgebraic(P = 10, V = 120)
expected_vsi_minimal = VSIMinimal(τ_P = 1, τ_Q = 2, K_P = 3, K_Q = 4, V_r = 5, P = 6, Q = 7)
expected_vsi_voltage_pt1 =
    VSIVoltagePT1(τ_v = 1, τ_P = 2, τ_Q = 3, K_P = 4, K_Q = 5, V_r = 6, P = 7, Q = 8)
expected_csi_minimal = CSIMinimal(I_r = 1)
expected_exponential_recov_load = ExponentialRecoveryLoad(
    P0 = 1,
    Q0 = 2,
    Nps = 3,
    Npt = 4,
    Nqs = 5,
    Nqt = 6,
    Tp = 7,
    Tq = 8,
    V0 = 9,
)
expected_fourth_order_eq_avr = FourthOrderEqGovernorExciterAVR(
    H = 1,
    P = 2,
    D = 3,
    Ω = 4,
    T_d_dash = 5,
    T_q_dash = 6,
    X_q_dash = 7,
    X_d_dash = 8,
    X_d = 9,
    X_q = 10,
    T_e = 11,
    T_a = 12,
    T_f = 13,
    K_e = 14,
    K_a = 15,
    K_f = 16,
    V_ref = 17,
    R_d = 18,
    T_sv = 19,
    T_ch = 20,
)

expected_static_line_array = StaticLine(from = 1, to = 2, Y = 0.1 + 5im)
expected_pimodel_line_array =
    PiModelLine(from = 1, to = 3, y = 0.1 + 5im, y_shunt_km = 200, y_shunt_mk = 300)

expected_static_line_dict = StaticLine(from = "bus1", to = "bus2", Y = 0.1 + 5im)
expected_pimodel_line_dict =
    PiModelLine(from = "bus1", to = "bus3", y = 0.1 + 5im, y_shunt_km = 200, y_shunt_mk = 300)

expected_nodes_array = [
    expected_swing_eq,
    expected_swing_eq_lvs,
    expected_fourth_order_eq,
    expected_slack,
    expected_pq,
    expected_pv,
    expected_vsi_minimal,
    expected_vsi_voltage_pt1,
    expected_csi_minimal,
    expected_exponential_recov_load,
    expected_fourth_order_eq_avr,
]
expected_lines_array = [expected_static_line_array, expected_pimodel_line_array]

expected_nodes_dict = OrderedDict(
    "bus1" => expected_swing_eq,
    "bus2" => expected_swing_eq_lvs,
    "bus3" => expected_fourth_order_eq,
    "bus4" => expected_slack,
    "bus5" => expected_pq,
    "bus6" => expected_pv,
    "bus7" => expected_vsi_minimal,
    "bus8" => expected_vsi_voltage_pt1,
    "bus9" => expected_csi_minimal,
    "bus10" => expected_exponential_recov_load,
    "bus11" => expected_fourth_order_eq_avr,
)
expected_lines_dict = OrderedDict("branch1" => expected_static_line_dict, "branch2" => expected_pimodel_line_dict)

array_grid = read_powergrid(joinpath(@__DIR__, "array_grid.json"), Json)
@test array_grid.nodes == expected_nodes_array
@test array_grid.lines == expected_lines_array

dict_grid = read_powergrid(joinpath(@__DIR__, "dict_grid.json"), Json)
@test dict_grid.nodes == expected_nodes_dict
@test dict_grid.lines == expected_lines_dict

target_dir = "../target"
array_export_file = joinpath(target_dir, "array_grid_export.json")
dict_export_file = joinpath(target_dir, "dict_grid_export.json")
mkpath(target_dir)

write_powergrid(array_grid, array_export_file, Json)
write_powergrid(dict_grid, dict_export_file, Json)

# complete the cycle and read back in
array_grid = read_powergrid(array_export_file, Json)
dict_grid = read_powergrid(dict_export_file, Json)

@test array_grid.nodes == expected_nodes_array
@test array_grid.lines == expected_lines_array

@test dict_grid.nodes == expected_nodes_dict
@test dict_grid.lines == expected_lines_dict

@test typeof(dict_grid.nodes)<: OrderedDict 
@test typeof(dict_grid.lines)<: OrderedDict 

@test typeof(array_grid.nodes)<:Array
@test typeof(array_grid.lines)<:Array


@test_throws ArgumentError read_powergrid(
    joinpath(@__DIR__, "grid_with_invalid_type.json"),
    Json,
)
@test_throws UndefKeywordError read_powergrid(
    joinpath(@__DIR__, "grid_with_invalid_params.json"),
    Json,
)

# test parsing of parameter-arrays

P, Q, V = rand(3)
Bᵤ = 1im*rand(2)
Cᵤ, Gᵤ, Hᵤ = 1im*rand(3)
Bₓ = rand(2,2)
Cₓ = rand(2)
Gₓ = rand(2)
Hₓ = rand(2)

nodes = [
    NormalForm(P=P, Q=Q, V=V, Bᵤ=Bᵤ, Cᵤ=Cᵤ, Gᵤ=Gᵤ, Hᵤ=Hᵤ, Bₓ=Bₓ, Cₓ=Cₓ, Gₓ=Gₓ, Hₓ=Hₓ)
]

pg = PowerGrid(nodes,[])

file = joinpath(@__DIR__,"testgrid.json")
write_powergrid(pg, file, Json)
pg_read = read_powergrid(file, Json)
@test pg_read.nodes == pg.nodes
rm(file)