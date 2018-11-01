
include("testing_base.jl")

################################################################################
# define variables as SymPy symbols for all node tests
@syms u i
@syms t real=true
# derivaed dynamic variables
v = abs(u)
s = u*conj(i)
p, q = real(s), imag(s)
################################################################################

@testset "PQAlgebraic" begin
@syms S
pq_par = PQAlgebraic(S=S)
pq_dyn = construct_node_dynamics(pq_par)
@test pq_par === parametersof(pq_dyn)
dint = []; int = []
@test pq_dyn.ode_dynamics.rhs(dint, u, i, int, t) == S - s
@test internalsymbolsof(pq_dyn) == []
@test internaldsymbolsof(pq_dyn) == []

@syms du
ints = PowerDynBase.ODEVariable(val = int)
us = PowerDynBase.ODEVariable(val = [u])
pq_dyn(1, us, i, ints, t)
@test us.ddt[1] == S - s
end

@testset "PVAlgebraic" begin
@syms P real=true
@syms V positive=true
pv_dyn = construct_node_dynamics(PVAlgebraic(P=P, V=V))
dint = []; int = []
@test pv_dyn.ode_dynamics.rhs(dint, u, i, int, t) == (v-V) + im*(p-P)
@test internalsymbolsof(pv_dyn) == []
@test internaldsymbolsof(pv_dyn) == []
end

@testset "SlackAlgebraic" begin
@syms U
slack_dyn = construct_node_dynamics(SlackAlgebraic(U=U))
dint = []; int = []
@test slack_dyn.ode_dynamics.rhs(dint, u, i, int, t) == u - U
@test internalsymbolsof(slack_dyn) == []
@test internaldsymbolsof(slack_dyn) == []
end

@testset "Swing Equations" begin
@syms H D positive=true
@syms P Ω real=true
@syms omega domega real=true
swing_par = SwingEq(H=H, P=P, D=D, Ω=Ω)
swing_dyn = construct_node_dynamics(swing_par)
@test swing_par === parametersof(swing_dyn)
dint = [domega]; int = [omega]; int_test = copy(int)
@test swing_dyn.rhs(dint, u, i, int, t) == u*im*omega
@test expand.(dint) == expand.([(P - D*omega - p)*2PI*Ω/H])
@test internalsymbolsof(swing_dyn) == [:ω]
@test internaldsymbolsof(swing_dyn) == [:dω]

@syms du
ints = PowerDynBase.ODEVariable(val = int, ddt = dint)
us = PowerDynBase.ODEVariable(val = [u], ddt = [du])
swing_dyn(1, us, i, ints, t)
@test us.ddt[1] == u*im*omega
@test expand.(ints.ddt) == expand.([(P - D*omega - p)*2PI*Ω/H])

@syms V Γ positive=true
swing_lvs_dyn = construct_node_dynamics(SwingEqLVS(H=H, P=P, D=D, Ω=Ω, Γ=Γ, V=V))
dint = [domega]; int = [omega]; int_test = copy(int)
@test swing_lvs_dyn.rhs(dint, u, i, int, t) == u*im*omega - u/v * Γ * (v-V)
@test expand.(dint) == expand.([(P - D*omega - p)*2PI*Ω/H])
@test internalsymbolsof(swing_lvs_dyn) == [:ω]
@test internaldsymbolsof(swing_lvs_dyn) == [:dω]

algebraic_swing_dyn = convert(PowerDynBase.AlgebraicNodeDynamics, swing_dyn)
@test swing_par === parametersof(algebraic_swing_dyn)
@syms out_omega real=true
@syms du
int_out = [out_omega]; dint = [domega]; int = [omega]
@test algebraic_swing_dyn.root(int_out, du,  dint, u, i, int, t) == u*im*omega - du
@test expand.(int_out) == expand.([(P - D*omega - p)*2PI*Ω/H]) .- dint
@test internalsymbolsof(algebraic_swing_dyn) == [:ω]
@test internaldsymbolsof(algebraic_swing_dyn) == [:dω]
@test internaloutsymbolsof(algebraic_swing_dyn) == [:outω] # automatically generated symbol
end

@testset "FourthEq" begin
@syms H  D  Ω  T_d_dash T_q_dash X_q_dash X_d_dash X_d X_q positive=true
@syms P  E_f real=true
@syms omega domega theta dtheta real=true
fourth_dyn = construct_node_dynamics(FourthEq(H=H, P=P, D=D, Ω=Ω, E_f=E_f, T_d_dash=T_d_dash ,T_q_dash=T_q_dash ,X_q_dash=X_q_dash ,X_d_dash=X_d_dash,X_d=X_d, X_q=X_q))
dint = [dtheta,domega]; int = [theta,omega]; int_test = copy(int)
du = fourth_dyn.rhs(dint, u, i, int, t)
i_c = 1im*i*exp(-1im*theta)
e_c = 1im*u*exp(-1im*theta)
p = real(u * conj(i))
e_d = real(e_c)
e_q = imag(e_c)
i_d = real(i_c)
i_q = imag(i_c)
de_d = (1 / T_q_dash)* (- e_d + (X_q - X_q_dash)* i_q)
de_q = (1 / T_d_dash)* (- e_q - (X_d - X_d_dash) * i_d + E_f)# sign error?
de_c = de_d + 1im*de_q
@test expand(dint[1]) == omega
@test du == -1im*de_c*exp(1im*theta)+ u*1im*omega
@test expand(dint[2]) == expand((P - D*omega - p- (X_q_dash - X_d_dash)*i_d* i_q)*2PI*Ω/H)
end

@testset "VSIMinimal" begin
@syms τ_P τ_Q K_P K_Q positive=true
@syms P Q V_r real=true
@syms omega domega real=true
VSIMindyn = construct_node_dynamics(VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q))
dint = [domega]; int = [omega]; int_test = copy(int)
du = VSIMindyn.rhs(dint, u, i, int, t)
p = real(u * conj(i))
q = imag(u * conj(i))
v = abs(u)
dv = 1/τ_Q*(-v+ V_r- K_Q *(q-Q))
@test du == u * 1im * omega + dv*(u/v)
@test expand.(dint) == expand.([1/τ_P*(-omega-K_P*(p-P))])
end

@testset "VSIVoltagePT1" begin
@syms τ_v τ_P τ_Q K_P K_Q positive=true
@syms P Q V_r real=true
@syms q_m dq_m omega domega real=true
VSIdyn = construct_node_dynamics(VSIVoltagePT1(τ_v=τ_v,τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q))
dint = [domega,dq_m]; int = [omega,q_m]; int_test = copy(int)
du = VSIdyn.rhs(dint, u, i, int, t)
p = real(u * conj(i))
q = imag(u * conj(i))
v = abs(u)
dv = 1/τ_v*(-v+V_r - K_Q*(q_m-Q))
@test du == u * 1im * omega + dv*(u/v)
@test expand.(dint[1]) == expand.(1/τ_P*(-omega-K_P*(p-P)))
@test expand.(dint[2]) == expand.(1/τ_Q*(q-q_m))
end
