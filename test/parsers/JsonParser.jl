using PowerDynamics: read, write, Json
using Test: @test

expected_swing_eq = SwingEq(H=1, P=2, D=3, Ω=12)
expected_swing_eq_lvs = SwingEqLVS(H=1, P=2, D=3, Ω=12, Γ=13, V=14)
expected_fourth_order_eq = FourthOrderEq(H=1, P=2, D=3, Ω=12, E_f=13, T_d_dash=14 ,T_q_dash=15 ,X_q_dash=16 ,X_d_dash=17,X_d=18, X_q=19)
expected_slack = SlackAlgebraic(U=10)

expected_static_line = StaticLine(from=1, to=2, Y=0.1+5im)
expected_pimodel_line = PiModelLine(from=1, to=3, y=0.1+5im, y_shunt_km=200, y_shunt_mk=300)

expected_nodes  = [expected_swing_eq, expected_swing_eq_lvs, expected_fourth_order_eq, expected_slack]
expected_lines = [expected_static_line, expected_pimodel_line]

grid = read(joinpath(@__DIR__, "grid.json"), Json)
@test grid.nodes == expected_nodes
@test grid.lines == expected_lines

target_dir = "../target"
export_file = joinpath(target_dir,"grid_export.json")
mkpath(target_dir)

write(grid, export_file, Json)

# complete the cycle and read back in
grid = read(export_file, Json)

@test grid.nodes == expected_nodes
@test grid.lines == expected_lines
