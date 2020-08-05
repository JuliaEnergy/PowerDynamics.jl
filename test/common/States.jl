using Test: @test, @testset, @test_throws
using PowerDynamics: SlackAlgebraic, SwingEqLVS, PowerGrid, State, startindex, StateError,StaticLine,systemsize
using LightGraphs: SimpleGraph, add_edge!, edges


Y = 0 + 5*im
graph = SimpleGraph(3)
nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
add_edge!(graph, 1, 2)
add_edge!(graph, 1, 3)
lines = [StaticLine(from=e.src, to=e.dst, Y=Y) for e in edges(graph)]
grid = PowerGrid(graph, nodes, lines)
u_Sl,u_Sw1,u_Sw2 = rand(ComplexF64, 3)
omega1,omega2 = rand(1.0:10.0, 2)
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])

nodes_dict = Dict("bus1"=>SlackAlgebraic(U=1), "bus2"=> SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), "bus3"=> SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1))
lines_dict = Dict("line1"=>StaticLine(from="bus1", to="bus2", Y=Y),"line2"=>StaticLine(from="bus1", to="bus3",Y=Y))
grid_fromdict = PowerGrid(graph, nodes_dict, lines_dict)
state_fromdict = State(grid_fromdict, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])




@test startindex(nodes_dict, "bus1") == 0
@test startindex(nodes_dict, "bus2") == 2
@test startindex(nodes_dict, "bus3") == 5


@test startindex(nodes, 1) == 0
@test startindex(nodes, 2) == 2
@test startindex(nodes, 3) == 5

@test state[1, :u] == complex(u_Sl)
@test state[2, :u] == complex(u_Sw1)
@test state[1:2, :u] == [complex(u_Sl), complex(u_Sw1)]

@test state_fromdict["bus1", :u] == complex(u_Sl)
@test state_fromdict["bus2", :u] == complex(u_Sw1)
@test state_fromdict[["bus1","bus2"], :u] == [complex(u_Sl), complex(u_Sw1)]

@test state_fromdict[:, :u] == [complex(u_Sl), complex(u_Sw1),complex(u_Sw2)]
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

@test state_fromdict["bus2", :v] == complex(u_Sw1) |> abs
@test state_fromdict["bus3", :v] == complex(u_Sw2) |> abs
@test state_fromdict["bus2", :var, 3] == omega1
@test state_fromdict["bus2", :ω] == omega1
@test state_fromdict["bus3", :var, 3] == omega2
@test state_fromdict["bus3", :ω] == omega2
@test state_fromdict[["bus2","bus3"], :var, 3] == [omega1, omega2]
@test state_fromdict[["bus2","bus3"], :ω] == [omega1, omega2]
#@test state_fromdict[["bus2","bus3"], :var, [3, 3]] == [omega1, omega2]


@test_throws BoundsError state[length(grid.nodes)+1, :u]
@test_throws BoundsError state[2, :var, 4]
@test_throws BoundsError state[1, :var, 3]
@test_throws BoundsError state[length(grid.nodes)+1, :var, 3]
@test_throws StateError state[1, :ω]

@test_throws BoundsError state_fromdict[length(grid.nodes)+1, :u]
@test_throws BoundsError state_fromdict["bus2", :var, 4]
@test_throws BoundsError state_fromdict["bus1", :var, 3]
@test_throws StateError state_fromdict["bus1", :ω]

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

# exchange omega
state_fromdict["bus2", :var, 3] = omega2
@test state_fromdict[["bus2","bus3"], :var, 3] == [omega2, omega2]
state_fromdict[["bus2","bus3"], :var, 3] = [omega1, omega2]
@test state_fromdict[["bus2","bus3"], :var, 3] == [omega1, omega2]
state_fromdict["bus2", :ω] = omega2
@test state_fromdict[["bus2","bus3"], :ω] == [omega2, omega2]
state_fromdict[["bus2","bus3"], :ω] = [omega1, omega2]
@test state_fromdict[["bus2","bus3"], :ω] == [omega1, omega2]

# modify the v as test
v = rand(1:10)
state[1, :u] = u_Sw1
state[:, :v] = v
@test state[1, :v] ≈ v

# modify the v as test
v = rand(1:10)
state_fromdict["bus1", :u] = u_Sw1
state_fromdict[:, :v] = v
@test state_fromdict["bus1", :v] ≈ v

state_fromdict["bus1", :u] = u_Sw1
@test state_fromdict["bus1", :u] == complex(u_Sw1)
state_fromdict["bus2", :u] = u_Sl
@test state_fromdict["bus2", :u] == complex(u_Sl)
# change back
state_fromdict[:, :u] = [u_Sl, u_Sw1, u_Sw2]
@test state_fromdict[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]


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

state_from_dict = State(grid_fromdict, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), omega1, real(u_Sw2), imag(u_Sw2), omega2])
@test state_from_dict["bus1", :φ] ≈ φ_Sl
@test state_from_dict["bus2", :φ] ≈ φ_Sw1
@test state_from_dict["bus3", :φ] ≈ φ_Sw2

i_1 = -Y * (u_Sw1 - u_Sl) -Y * (u_Sw2 - u_Sl)
i_2 = Y * (u_Sw1 - u_Sl)
i_3 = Y * (u_Sw2 - u_Sl)
@test state[1, :i] == i_1
@test state[2, :i] == i_2

@test state_from_dict["bus1", :i] == i_1
@test state_from_dict["bus2", :i] == i_2

@test state[:, :i] == [i_1, i_2, i_3]
@test state[2, :iabs] == abs(i_2)
s_2 =  state[2, :u] * conj(i_2)
@test state[2, :s] == s_2
@test state[2, :p] ==  real(s_2)
@test state[2, :q] ==  imag(s_2)
@test state[:, :s] ==  state[:, :u] .* conj.([i_1, i_2, i_3])
@test state[:, :p] ==  real(state[:, :u] .* conj.([i_1, i_2, i_3]))
@test state[:, :q] ==  imag(state[:, :u] .* conj.([i_1, i_2, i_3]))

@test state_from_dict[:, :i] == [i_1, i_2, i_3]
@test state_from_dict["bus2", :iabs] == abs(i_2)
s_2 =  state_from_dict["bus2", :u] * conj(i_2)
@test state_from_dict["bus2", :s] == s_2
@test state_from_dict["bus2", :p] ==  real(s_2)
@test state_from_dict["bus2", :q] ==  imag(s_2)
@test state_from_dict[:, :s] ==  state[:, :u] .* conj.([i_1, i_2, i_3])
@test state_from_dict[:, :p] ==  real(state[:, :u] .* conj.([i_1, i_2, i_3]))
@test state_from_dict[:, :q] ==  imag(state[:, :u] .* conj.([i_1, i_2, i_3]))

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

s1_fromdict = State(grid_fromdict, ones(systemsize(grid)))
s2_fromdict = State(grid_fromdict, ones(systemsize(grid)))
@test s1_fromdict == s2_fromdict
@test s1_fromdict ≈ s2_fromdict
@test ((s1_fromdict + s2_fromdict).vec .== 2) |> all
@test ((s1_fromdict - s2_fromdict).vec .== 0) |> all
@test ((-s1_fromdict).vec .== -1) |> all
@test ((2*s1_fromdict).vec .== 2) |> all
@test ((s1_fromdict*2).vec .== 2) |> all
@test ((s1_fromdict/2).vec .≈ 0.5) |> all

@test convert(Vector, s1) == ones(Float64, 8)
@test convert(Array, s1) == ones(Float64, 8)

@test convert(Vector, s1_fromdict) == ones(Float64, 8)
@test convert(Array, s1_fromdict) == ones(Float64, 8)
