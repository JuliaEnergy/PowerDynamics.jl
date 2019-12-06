using PowerDynamics
import Base: @__doc__ #-> LoadError: UndefVarError: @__doc__ not defined
import PowerDynamics: AbstractNode #-> UndefVarError: AbstractNode not defined


@DynamicNode MySwingEq(P, H, D, Ω) begin
PowerDynamics.MassMatrix(;m_u = true, m_int = [1, 0])
end begin
@assert D > 0 "damping (D) should be >0"
@assert H > 0 "inertia (H) should be >0"
Ω_H = Ω * 2pi / H
end [[ω, dω]] begin
p = real(u * conj(i_c))
dϕ = ω # dϕ is only a temp variable that Julia should optimize out
du = u * im * dϕ
dω = (P - D*ω - p)*Ω_H
end

showdefinition(stdout, MySwingEq) |> println

swing_node = MySwingEq(P=1., H=1., D=0.1, Ω=50.)
slack = SlackAlgebraic(U=1.)

nodes = [swing_node, slack]

lines = [StaticLine(from=1, to=2, Y=1. + 5. * im)]

pg = PowerGrid(nodes, lines)

# The methods created by the macro are working here...

dimension(swing_node)
symbolsof(swing_node)

# The following two lines both fail with:
#   MethodError: no method matching dimension(::MySwingEq)
# This leads to subsequent errors e.g. in find_operationpoint or find_valid_initial_condition.

systemsize(pg)
