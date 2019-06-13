
include("testing_base.jl")

nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
numnodes = length(nodes)
LY = SymbolMatrix("LY", numnodes)
grid = GridDynamics(nodes, LY, skip_LY_check=true)
@syms u_Sl u_Sw1 u_Sw2
@syms omega1 omega2 real=true
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), real(u_Sw2), imag(u_Sw2), omega1, omega2])

@test state[1, :u] == complex(u_Sl)
@test state[2, :u] == complex(u_Sw1)
@test state[1:2, :u] == [complex(u_Sl), complex(u_Sw1)]
@test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]
@test state[2, :v] == complex(u_Sw1) |> abs
@test state[3, :v] == complex(u_Sw2) |> abs
@test state[2, :int, 1] == omega1
@test state[2, :ω] == omega1
@test state[3, :int, 1] == omega2
@test state[3, :ω] == omega2
@test state[2:3, :int, 1] == [omega1, omega2]
@test state[2:3, :ω] == [omega1, omega2]
@test state[2:3, :int, [1, 1]] == [omega1, omega2]
@test state[1, :i] == (LY * complex.([u_Sl, u_Sw1, u_Sw2]))[1]
@test state[:, :i] == LY * complex.([u_Sl, u_Sw1, u_Sw2])
@test state[1, :iabs] == abs((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[1])
@test state[:, :iabs] == abs.((LY * complex.([u_Sl, u_Sw1, u_Sw2])))
@test state[2, :s] ==  state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2])
@test state[2, :p] ==  real(state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2]))
@test state[2, :q] ==  imag(state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2]))
@test state[:, :s] ==  state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2])))
@test state[:, :p] ==  real(state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2]))))
@test state[:, :q] ==  imag(state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2]))))

@test_throws BoundsError state[numnodes+1, :u]
@test_throws BoundsError state[2, :int, 2]
@test_throws BoundsError state[1, :int, 1]
@test_throws BoundsError state[numnodes+1, :int, 1]
@test_throws StateError state[1, :ω]

# exchange the u as test
state[1, :u] = u_Sw1
@test state[1, :u] == complex(u_Sw1)
state[2, :u] = u_Sl
@test state[2, :u] == complex(u_Sl)
# change back
state[:, :u] = [u_Sl, u_Sw1, u_Sw2]
@test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]
@test_throws BoundsError state[numnodes+1, :u] = u_Sw1
# exhange omega
state[2, :int, 1] = omega2
@test state[2:3, :int, 1] == [omega2, omega2]
state[2:3, :int, 1] = [omega1, omega2]
@test state[2:3, :int, 1] == [omega1, omega2]
state[2, :ω] = omega2
@test state[2:3, :ω] == [omega2, omega2]
state[2:3, :ω] = [omega1, omega2]
@test state[2:3, :ω] == [omega1, omega2]
@test_throws BoundsError state[1, :int, 1] = omega2
@test_throws StateError state[1, :ω] = omega2

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
state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), real(u_Sw2), imag(u_Sw2), omega1, omega2])
@test state[1, :φ] ≈ φ_Sl
@test state[2, :φ] ≈ φ_Sw1
@test state[3, :φ] ≈ φ_Sw2

## binary operations ##
s1 = State(grid, ones(SystemSize(grid)))
s2 = State(grid, ones(SystemSize(grid)))
@test s1 == s1
@test s1 ≈ s2
@test (PowerDynBase.BaseState(s1 + s2).vec .== 2) |> all
@test (PowerDynBase.BaseState(s1 - s2).vec .== 0) |> all
@test (PowerDynBase.BaseState(-s1).vec .== -1) |> all
@test (PowerDynBase.BaseState(2*s1).vec .== 2) |> all
@test (PowerDynBase.BaseState(s1*2).vec .== 2) |> all
@test (PowerDynBase.BaseState(s1/2).vec .≈ 0.5) |> all

@test convert(PowerDynBase.BaseState, s1) === s1.base
@test convert(Vector, s1) == ones(Float64, 8)
@test convert(Array, s1) == ones(Float64, 8)
