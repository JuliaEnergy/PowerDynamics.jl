using Test: @test, @testset
using PowerDynBase: SlackAlgebraic, SwingEqLVS, PowerGrid, State
using SymPy: @syms, simplify
using LightGraphs: SimpleGraph, add_edge!, edges

Y = 0 + 5*im
graph = SimpleGraph(3)
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
add_edge!(graph, 1, 2)
add_edge!(graph, 1, 3)
lines = [StaticLine(Y=Y) for e in edges(graph)]
grid = PowerGrid(graph, nodes, lines)
@syms u_Sl u_Sw1 u_Sw2
@syms omega1 omega2 real=true
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])

@test state[1, :u] == complex(u_Sl)
@test state[2, :u] == complex(u_Sw1)
@test state[1:2, :u] == [complex(u_Sl), complex(u_Sw1)]

@test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]

@test state[2, :v] == complex(u_Sw1) |> abs
@test state[3, :v] == complex(u_Sw2) |> abs
@test state[2, :var, 3] == omega1
@test state[2, :ω] == omega1
@test state[3, :var, 3] == omega2
@test state[3, :ω] == omega2
@test state[2:3, :var, 3] == [omega1, omega2]
@test state[2:3, :ω] == [omega1, omega2]
@test state[2:3, :var, [3, 3]] == [omega1, omega2]

@test_throws BoundsError state[length(grid.nodes)+1, :u]
@test_throws BoundsError state[2, :var, 4]
@test_throws BoundsError state[1, :var, 3]
@test_throws BoundsError state[length(grid.nodes)+1, :var, 3]
@test_throws StateError state[1, :ω]

# exchange the u as test
state[1, :u] = u_Sw1
@test state[1, :u] == complex(u_Sw1)
state[2, :u] = u_Sl
@test state[2, :u] == complex(u_Sl)
# change back
state[:, :u] = [u_Sl, u_Sw1, u_Sw2]
@test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]
@test_throws BoundsError state[length(grid.nodes)+1, :u] = u_Sw1

# exchange omega
state[2, :var, 3] = omega2
@test state[2:3, :var, 3] == [omega2, omega2]
state[2:3, :var, 3] = [omega1, omega2]
@test state[2:3, :var, 3] == [omega1, omega2]
state[2, :ω] = omega2
@test state[2:3, :ω] == [omega2, omega2]
state[2:3, :ω] = [omega1, omega2]
@test state[2:3, :ω] == [omega1, omega2]

# modify the v as test
@syms v positive=true
state[1, :u] = u_Sw1
state[:, :v] = v
@test state[1, :v] |> simplify == v

## getindex for angle (numerically) ##
v_Sl = rand()
v_Sw1 = rand()
v_Sw2 = rand()
φ_Sl = rand()
φ_Sw1 = rand()
φ_Sw2 = rand()
u_Sl = v_Sl*exp(im*φ_Sl)
u_Sw1 = v_Sw1*exp(im*φ_Sw1)
u_Sw2 = v_Sw2*exp(im*φ_Sw2)
omega1 = 3.
omega2 = 5.
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])
@test state[1, :φ] ≈ φ_Sl
@test state[2, :φ] ≈ φ_Sw1
@test state[3, :φ] ≈ φ_Sw2

i_1 = -Y * (u_Sw1 - u_Sl) -Y * (u_Sw2 - u_Sl)
i_2 = Y * (u_Sw1 - u_Sl)
i_3 = Y * (u_Sw2 - u_Sl)
@test state[1, :i] == i_1
@test state[2, :i] == i_2

@test state[:, :i] == [i_1, i_2, i_3]
@test state[2, :iabs] == abs(i_2)
s_2 =  state[2, :u] * conj(i_2)
@test state[2, :s] == s_2
@test state[2, :p] ==  real(s_2)
@test state[2, :q] ==  imag(s_2)
@test state[:, :s] ==  state[:, :u] .* conj.([i_1, i_2, i_3])
@test state[:, :p] ==  real(state[:, :u] .* conj.([i_1, i_2, i_3]))
@test state[:, :q] ==  imag(state[:, :u] .* conj.([i_1, i_2, i_3]))

## binary operations ##
s1 = State(grid, ones(systemsize(grid)))
s2 = State(grid, ones(systemsize(grid)))
@test s1 == s1
@test s1 ≈ s2
@test ((s1 + s2).vec .== 2) |> all
@test ((s1 - s2).vec .== 0) |> all
@test ((-s1).vec .== -1) |> all
@test ((2*s1).vec .== 2) |> all
@test ((s1*2).vec .== 2) |> all
@test ((s1/2).vec .≈ 0.5) |> all

@test convert(Vector, s1) == ones(Float64, 8)
@test convert(Array, s1) == ones(Float64, 8)
