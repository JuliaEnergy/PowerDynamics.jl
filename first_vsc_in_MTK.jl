# First attempt at VSC in BlockSystems (nee IOSystems)

using PowerDynamics: SlackAlgebraic,GFI,PiModelLine,PowerGrid,find_operationpoint
using OrderedCollections: OrderedDict

##
τ_P, τ_Q, K_P, K_Q = rand(1:10, 4)
P,Q,V_r = rand(1.:10.,3)

buses=OrderedDict(
    "bus1"=> GFI(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_refP=V_r,V_refQ=V_r,P_ref=P,Q_ref=Q),
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

##
powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,5.)

# simulating a frequency perturbation at node 1
#fault1 = ChangeInitialConditions(node="bus1", var=:ω, f=Inc(0.2))
#solution1 = simulate(fault1, powergrid, operationpoint, timespan)
