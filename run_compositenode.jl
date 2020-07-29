using PowerDynamics: CompositeNode, CSIMinimal, VSIMinimal, PQAlgebraic, construct_vertex, symbolsof, SlackAlgebraic,PiModelLine,PowerGrid,find_operationpoint
using PowerDynamics: PowerPerturbation,PQAlgebraic
# constant voltage
#VSI = VSIMinimal(τ_P=1.,τ_Q=1.,K_P=1.,K_Q=1.,V_r=1.,P=1.,Q=1.)

τ_P, τ_Q, K_P, K_Q = rand_real(4)
P,Q,V_r =rand_real(3)
VSI = VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q)

# constant current
I_r =rand(1.0:10.0)+1im*rand(1.0:10.0)
CSI = CSIMinimal(I_r=I_r)
CSI2 = CSIMinimal(I_r=I_r)

# constant power
P,Q=rand_real(2)
PQ = PQAlgebraic(P=P,Q=Q)
PQ2 = PQAlgebraic(P=P,Q=Q)

# Test vertex with voltage dynamic
comp_node = CompositeNode(CurrentNodes=[CSI, CSI2], PowerNodes=[PQ, PQ2], VoltageNode=VSI)
comp_vertex = construct_vertex(comp_node)

node_list=[
    SlackAlgebraic(U=1),
    PQAlgebraic(P=1,Q=1),
    comp_node]

line_list=[
    PiModelLine(from=1, to=2, y=4.99913-15.2631im, y_shunt_km=1/2*0.0528*1im, y_shunt_mk=1/2*0.0528*1im),
    PiModelLine(from=2, to=3, y=4.99913-15.2631im, y_shunt_km=1/2*0.0528*1im, y_shunt_mk=1/2*0.0528*1im)]

final_time=5.
timespan= (0.,final_time);
disturbed_node=2;
powergrid = PowerGrid(node_list, line_list);
operationpoint = find_operationpoint(powergrid);
pd = PowerPerturbation(node_number=disturbed_node,fraction=0.5, tspan_fault=(0.5,5));
#include("plotting2.jl");

sol_pd,result_pd = simulate(pd,powergrid,operationpoint,timespan)
