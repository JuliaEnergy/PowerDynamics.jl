# Now let's test whether this works in PD...
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

buses=OrderedDict(
    "bus1"=> gfi,
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,0.1)

# simulating a voltage perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)
nothing #hide

#=
it works 🎉
## Benchmark
We can also run a quick benchmark of both node functions:
=#
para = Dict(:τ_v=>rand(),:τ_P=>rand(), :τ_Q=>rand(),
            :K_P=>rand(), :K_Q=>rand(), :V_r=>rand(),
            :P=>rand(), :Q=>rand(), :ω_r=>0.0)

node_bs = GFI(;para...)
f_bs = construct_vertex(node_bs).f!
para_pd = delete!(copy(para), :ω_r)
node_pd = VSIVoltagePT1(; para_pd...)
f_pd = construct_vertex(node_pd).f!

edges = [randn(4) for i in 1:4]
t = rand()

## chose random initial state and account for initial ω in PD node
x_bs = randn(4)
x_pd = copy(x_bs)
x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])
dx = similar(x_bs)

using BenchmarkTools
@btime $f_bs($dx, $x_bs, $edges, nothing, $t)
@btime $f_pd($dx, $x_pd, $edges, nothing, $t)
# it seems like the `IONode` is even a bit faster.
