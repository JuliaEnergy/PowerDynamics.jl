using Pkg;
Pkg.instantiate();
cd(@__DIR__);
using PowerDynamics;
using OrdinaryDiffEq: ODEProblem, Rodas4, solve;
include("plotting.jl");
#include("helpers.jl");

busses = [VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  SwingEqLVS(H=1, P=1, D=1, Ω=2π*50,Γ=0.1,V=1.),
  #VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  VSIMinimal(;τ_P=1.0,τ_Q=0.0001,K_P=0.2,K_Q=0.002,V_r=1.0,P=1.0,Q=0.5),
  PQAlgebraic(P=-1.0,Q=-0.5),
  PQAlgebraic(P=-1.0,Q=-0.5),
  PQAlgebraic(P=-1.0,Q=-0.5)];
  #RLC_Load(R,L,C),RLC_Load(R,L,C),RLC_Load(R,L,C)]


lines = [
  PiModelLine(;from=1, to=4, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=5, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=3, to=6, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=2, to=3, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.),
  PiModelLine(;from=1, to=2, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.)];


pg = PowerGrid(busses, lines);
nsc = NodeShortCircuit(node = 2, R = 0.1, sc_timespan=(1.,2.));
ic_guess = ones(systemsize(pg)) + 0. * randn(systemsize(pg))
ic_guess = find_valid_initial_condition(pg,ic_guess)
problem = ODEProblem(rhs(pg),ic_guess,(0.,200.))
sol = solve(problem, Rodas4(autodiff=false))
op_point = sol[end];

sol_nc = simulate(nsc,pg,op_point,(0.,100.))
plot_res(sol_nc,pg,2)
