
include("testing_base.jl")

struct DummyNodeDynamics{N <: PowerDynBase.AbstractNodeParameters} <: PowerDynBase.AbstractNodeDynamics{N} end
let
println("OUTPUT TESTS:")
@test_nowarn println(PQAlgebraic)
@test_nowarn println(SwingEq)
@test_nowarn println(PQAlgebraic(S=3+4im))
@test_nowarn println(SwingEq(H=1, P=2, D=3, Ω=4))
Base.showerror(stdout, NodeDynamicsError("My message."))
println()
Base.showerror(stdout, GridDynamicsError("My message."))
println()
LY = [1 -1; 0 0]
@test_throws GridDynamicsError PowerDynBase.checkLY(LY)
@test_throws GridDynamicsError PowerDynBase.checkLY(transpose(LY))

LY = [1 -1; -1 1] # now a valid admittance laplacian
@test_throws GridDynamicsError GridDynamics([DummyNodeDynamics{PowerDynBase.SwingEq}()], LY)

@test_throws UndefKeywordError SlackAlgebraic()
end
