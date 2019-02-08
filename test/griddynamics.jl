
include("testing_base.jl")

@testset "OrdinaryGridDynamics" begin
@syms t real=true
num = 2
pars = [SwingEq(H=symbols("H_$i", positive=true), P=symbols("P_$i", real=true), D=symbols("D_$i", positive=true), Ω=symbols("Omega_$i", real=true)) for i=1:num ]
dyns = construct_node_dynamics.(pars)
# LY = [symbols("LY_$i$j") for i in 1:num, j in 1:num]
LY = SymbolMatrix("LY", num)
grid_dyn = GridDynamics(pars, LY, skip_LY_check=true)
@test grid_dyn isa PowerDynBase.OrdinaryGridDynamics
@test pars == grid_dyn |> Nodes .|> parametersof

# evalute the grid rhs symbolically
us = [symbols("u_$i") for i in 1:num ]
omegas = [symbols("omega_$i",real=true) for i=1:num ]
xs = [splitcomplexsymbols(us); omegas]
dxs = [
        [symbols("d$u")[1] for u in us]|> splitcomplexsymbols;
        [symbols("d$omega") for omega in omegas]
        ]
dxs_test = copy(dxs)
grid_dyn(dxs, xs, nothing, t)
dus = dxs[1:2*num] |> mergecomplexsymbols
domegas = dxs[2*num+1:end]

# Rebuild the grid rhs with the assumption that the node dynamics is working
# already (as it was tested separately before).
us = [symbols("u_$i") for i in 1:num ] .|> complex
us_var = PowerDynBase.ODEVariable(us)
omegas = [symbols("omega_$i",real=true) for i=1:num ]
omegas_var = PowerDynBase.ODEVariable(omegas)
t_i = LY*us
map(i -> dyns[i](i, us_var, t_i, view(omegas_var, i:i), nothing), 1:num)
@test all(complex.(dus) .== complex.(us_var.ddt)) # check whether the voltage differentials match
@test all(complex.(domegas) .== complex.(omegas_var.ddt)) # check whether the internal differentials match
end


@testset "OrdinaryGridDynamicsWithMass" begin
@syms t real=true
swing_num = 2
pq_num = 2
num = pq_num + swing_num
pars = [
    [SwingEq(H=symbols("H_$i", positive=true), P=symbols("P_$i", real=true), D=symbols("D_$i", positive=true), Ω=symbols("Omega_$i", real=true)) for i=1:swing_num ];
    [PQAlgebraic(S=symbols("S_$i")) for i in 1:pq_num]
    ]
dyns = construct_node_dynamics.(pars)
LY = SymbolMatrix("LY", num)
grid_dyn = GridDynamics(pars, LY, skip_LY_check=true)
@test typeof(grid_dyn) === PowerDynBase.OrdinaryGridDynamicsWithMass
@test PowerDynBase.masses(grid_dyn) == [true, true, true, true, false, false, false, false, true, true]
@test pars == grid_dyn |> Nodes .|> parametersof

# evalute the grid rhs symbolically
dyns = convert(Array{OrdinaryNodeDynamicsWithMass}, dyns)
us = [symbols("u_$i") for i in 1:num ]
omegas = [symbols("omega_$i",real=true) for i=1:swing_num ]
xs = [splitcomplexsymbols(us); omegas]
dxs = [
[symbols("d$u")[1] for u in us]|> splitcomplexsymbols;
[symbols("d$omega") for omega in omegas]
]
dxs_test = copy(dxs)
grid_dyn(dxs, xs, nothing, t)
dus = dxs[1:2*num] |> mergecomplexsymbols
domegas = dxs[2*num+1:end]

# Rebuild the grid rhs with the assumption that the node dynamics is working
# already (as it was tested separately before).
us = [symbols("u_$i") for i in 1:num ] .|> complex
us_var = PowerDynBase.ODEVariable(us)
omegas = [symbols("omega_$i",real=true) for i=1:swing_num ]
omegas_var = PowerDynBase.ODEVariable(omegas)
t_i = LY*us
map(i -> dyns[i](i, us_var, t_i, view(omegas_var, i:i), nothing), 1:swing_num);
map(i -> dyns[i](i, us_var, t_i, view(omegas_var,1:0), nothing), swing_num + 1:num);
@test all(complex.(dus) .== complex.(us_var.ddt)) # check whether the voltage differentials match
@test all(complex.(domegas) .== complex.(omegas_var.ddt)) # check whether the internal differentials match
end

@testset "OrdinaryGridDynamicsWithMass (numerically)" begin
nodes = [SwingEq(H=1, P=1, D=1, Ω=50), SwingEq(H=1, P=-1, D=1, Ω=50)]
LY = [im -im; -im im]
grid = GridDynamics(nodes, LY)
x = rand(SystemSize(grid))
dx = similar(x)
grid(dx, x, nothing, 0)
@test true # if the code runs until here, everything succeeded
end
