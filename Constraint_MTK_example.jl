using Revise
using PowerDynamics
using OrderedCollections: OrderedDict

gfi = GFI(τ_v = 0.005,    # time constant voltage control delay
          τ_P = 0.5,      # time constant active power measurement
          τ_Q = 0.5,      # time constant reactive power measurement
          K_P = 0.396,     # droop constant frequency droop
          K_Q = 0.198,     # droop constant voltage droop
          V_r = 1.0,   # reference/ desired voltage
          P   = 0.303, # active (real) power infeed
          Q   = 0.126, # reactive (imag) power infeed                .
          ω_r = 0.0)   # refrence/ desired frequency

pql = PQ_Load(P = -0.5, # active power set point
              Q = 0.0) # reactive power set point

buses=OrderedDict(
    "bus1"=> gfi,
    "bus2"=> SlackAlgebraic(U=1),
    "bus3"=> pql) # PQAlgebraic(P=0.5,Q=0.0))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2),
    "branch2"=> PiModelLine(from= "bus1", to = "bus3",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,0.1)

# simulating a voltage perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)
nothing #hide