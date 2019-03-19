
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

@testset "ExponetialRecoveryLoad" begin
@syms V0 Nps Npt Nqs Nqt Tp Tq positive=true
@syms P0 Q0 Pd Qd real=true
@syms x_p dx_p x_q dx_q real=true
ExpRec = construct_node_dynamics(ExponentialRecoveryLoad(P0=P0, Q0=Q0, Nps=Nps, Npt=Npt, Nqs=Nqs, Nqt=Nqt, Tp=Tp, Tq=Tq, V0=V0))
dint = [dx_p,dx_q]; int = [x_p,x_q]; int_test = copy(int)

du = ExpRec.ode_dynamics.rhs(dint, u, i, int, t)

@test expand(dint[1]) == expand((1/Tp) * (-x_p + P0*(v/V0)^Nps - P0*(v/V0)^Npt))
@test expand(dint[2]) == expand((1/Tq) * (-x_q + Q0*(v/V0)^Nqs - Q0*(v/V0)^Nqt))
@test expand(du) == expand((-p + x_p + P0*((v/V0)^Npt)) + im*(-q + x_q + Q0*((v/V0)^Nqt)))
end

@testset "FourthOrderEqGovernorExciterAVR" begin
@syms H  D  Ω  T_d_dash T_q_dash X_q_dash X_d_dash X_d X_q T_e T_a T_f K_a K_f V_ref R_d T_sv T_ch positive=true
@syms P K_e real=true
@syms omega domega theta dtheta e_f de_f v_r dv_r r_f dr_f P_sv dP_sv P_m dP_m real=true
fourth_dyn = construct_node_dynamics(FourthOrderEqGovernorExciterAVR(H=H, P=P, D=D, Ω=Ω, T_d_dash=T_d_dash ,T_q_dash=T_q_dash ,X_q_dash=X_q_dash ,X_d_dash=X_d_dash,X_d=X_d, X_q=X_q, T_e=T_e, T_a=T_a, T_f=T_f, K_e=K_e, K_a=K_a, K_f=K_f, V_ref=V_ref, R_d=R_d, T_sv=T_sv, T_ch=T_ch))
dint = [dtheta,domega,de_f,dv_r,dr_f,dP_sv,dP_m]; int = [theta,omega,e_f,v_r,r_f,P_sv,P_m]; int_test = copy(int)
du = fourth_dyn.rhs(dint, u, i, int, t)
i_c = 1im*i*exp(-1im*theta)
e_c = 1im*u*exp(-1im*theta)
p = real(u * conj(i))
e_d = real(e_c)
e_q = imag(e_c)
i_d = real(i_c)
i_q = imag(i_c)
V_mes = e_c - 1im*X_d_dash*i_c
de_d = (1 / T_q_dash)* (- e_d + (X_q - X_q_dash)* i_q)
de_q = (1 / T_d_dash)* (- e_q - (X_d - X_d_dash) * i_d + e_f)# sign error?
de_c = de_d + 1im*de_q
@test expand(dint[1]) == omega
@test du == -1im*de_c*exp(1im*theta)+ u*1im*omega
@test expand(dint[2]) == expand((P_m - D*omega - p- (X_q_dash - X_d_dash)*i_d* i_q)*2PI*Ω/H)
@test expand(dint[3]) == expand((1 / T_e) * ((- (K_e + (0.098*exp(0.55*e_f))) * e_f) + v_r))
@test expand(dint[4]) == expand((1 / T_a) * (- v_r + (K_a * r_f) - ((K_a * K_f)/T_f)*e_f + K_a*(V_ref - abs(V_mes))))
@test expand(dint[5]) == expand((1 / T_f) * (- r_f + ((K_f/T_f) * e_f)))
@test expand(dint[6]) == expand((1 / T_sv) * (-P_sv + P - (1/R_d)*(((omega+(Ω*2PI))/(Ω*2PI))-1)))
@test expand(dint[7]) == expand((1 / T_ch) * (-P_m  + P_sv))
end
