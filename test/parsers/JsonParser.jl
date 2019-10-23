using PowerDynamics: read_powergrid, write_powergrid, Json
using Test: @test

expected_swing_eq = SwingEq(H=1, P=2, D=3, Ω=12)
expected_swing_eq_lvs = SwingEqLVS(H=1, P=2, D=3, Ω=12, Γ=13, V=14)
expected_fourth_order_eq = FourthOrderEq(H=1, P=2, D=3, Ω=12, E_f=13, T_d_dash=14 ,T_q_dash=15 ,X_q_dash=16 ,X_d_dash=17,X_d=18, X_q=19)
expected_slack = SlackAlgebraic(U=10)
expected_pq = PQAlgebraic(P=0,Q=10)
expected_pv = PVAlgebraic(P=10, V=120)
expected_vsi_minimal = VSIMinimal(τ_P=1,τ_Q=2,K_P=3,K_Q=4,V_r=5,P=6,Q=7)
expected_vsi_voltage_pt1 = VSIVoltagePT1(τ_v=1,τ_P=2,τ_Q=3,K_P=4,K_Q=5,V_r=6,P=7,Q=8)
expected_csi_minimal = CSIMinimal(I_r=1)
expected_exponential_recov_load = ExponentialRecoveryLoad(P0=1, Q0=2, Nps=3, Npt=4, Nqs=5, Nqt=6, Tp=7, Tq=8, V0=9)
expected_fourth_order_eq_avr = FourthOrderEqGovernorExciterAVR(
H=1, P=2, D=3, Ω=4, T_d_dash=5 ,T_q_dash=6,
X_q_dash=7 ,X_d_dash=8,X_d=9, X_q=10, T_e=11, T_a=12, T_f=13, K_e=14, K_a=15,
K_f=16, V_ref=17, R_d=18, T_sv=19, T_ch=20)

expected_static_line = StaticLine(from=1, to=2, Y=0.1+5im)
expected_pimodel_line = PiModelLine(from=1, to=3, y=0.1+5im, y_shunt_km=200, y_shunt_mk=300)

expected_nodes  = [expected_swing_eq, expected_swing_eq_lvs, expected_fourth_order_eq, expected_slack,
expected_pq, expected_pv, expected_vsi_minimal, expected_vsi_voltage_pt1, expected_csi_minimal,
expected_exponential_recov_load, expected_fourth_order_eq_avr]
expected_lines = [expected_static_line, expected_pimodel_line]

grid = read_powergrid(joinpath(@__DIR__, "grid.json"), Json)
@test grid.nodes == expected_nodes
@test grid.lines == expected_lines

target_dir = "../target"
export_file = joinpath(target_dir,"grid_export.json")
mkpath(target_dir)

write_powergrid(grid, export_file, Json)

# complete the cycle and read back in
grid = read_powergrid(export_file, Json)

@test grid.nodes == expected_nodes
@test grid.lines == expected_lines


@test_throws ArgumentError read_powergrid(joinpath(@__DIR__, "grid_with_invalid_type.json"), Json)
@test_throws UndefKeywordError read_powergrid(joinpath(@__DIR__, "grid_with_invalid_params.json"), Json)
