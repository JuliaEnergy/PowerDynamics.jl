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

using Test
@testset "Compare IONode vs. PD Node" begin
    for i in 1:100
        para = Dict(
            :τ_v => rand(), # time constant voltage control delay
            :τ_P => rand(), # time constant active power measurement
            :τ_Q => rand(), # time constant reactive power measurement
            :K_P => rand(), # droop constant frequency droop
            :K_Q => rand(), # droop constant voltage droop
            :V_r => rand(), # reference/ desired voltage
            :P   => rand(), # active (real) power infeed
            :Q   => rand(), # reactive (imag) power infeed                .
            :ω_r => 0.0)    # refrence/ desired frequency

        ## create IONode
        node_bs = IONode(connected, para)
        f_bs = construct_vertex(node_bs).f!

        ## create PDNode
        ## the PD node does not have the explicit ω_r parameter
        para_pd = delete!(copy(para), :ω_r)
        nt = NamedTuple{Tuple(keys(para_pd))}(values(para_pd))
        node_pd = VSIVoltagePT1(; nt...)
        f_pd = construct_vertex(node_pd).f!

        ## create fake "edge data", 4 incoming, 4 outgooing with 4 states each
        es = [randn(4) for i in 1:4]
        ed = [randn(4) for i in 1:4]

        ## select random time
        t = rand()

        ## chose random initial state and account for initial ω in PD node
        x_bs = randn(4)
        x_pd = copy(x_bs)
        x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])

        ## create arrays for the results
        dx_bs = similar(x_bs)
        dx_pd = similar(x_pd)

        ## call both functions
        f_bs(dx_bs, x_bs, es, ed, nothing, t)
        f_pd(dx_pd, x_pd, es, ed, nothing, t)

        ## compare results
        ## we have to correct variable 3 of bs implementation to match dω
        @test dx_bs[1] ≈ dx_pd[1]                # u_r
        @test dx_bs[2] ≈ dx_pd[2]                # u_i
        @test - para[:K_P] * dx_bs[3] ≈ dx_pd[3] # ω
        @test dx_bs[4] ≈ dx_pd[4]                # q_filtered
    end
end
nothing #hide

#=
it works 🎉
## Benchmark
We can also run a quick benchmark of both node functions:
=#
para = Dict(:τ_v=>rand(),:τ_P=>rand(), :τ_Q=>rand(),
            :K_P=>rand(), :K_Q=>rand(), :V_r=>rand(),
            :P=>rand(), :Q=>rand(), :ω_r=>0.0)

node_bs = IONode(connected, para)
f_bs = construct_vertex(node_bs).f!
para_pd = delete!(copy(para), :ω_r)
nt = NamedTuple{Tuple(keys(para_pd))}(values(para_pd))
node_pd = VSIVoltagePT1(; nt...)
f_pd = construct_vertex(node_pd).f!

es = [randn(4) for i in 1:4]
ed = [randn(4) for i in 1:4]
t = rand()

## chose random initial state and account for initial ω in PD node
x_bs = randn(4)
x_pd = copy(x_bs)
x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])
dx = similar(x_bs)

using BenchmarkTools
@btime $f_bs($dx, $x_bs, $es, $ed, nothing, $t)
@btime $f_pd($dx, $x_pd, $es, $ed, nothing, $t)
# it seems like the `IONode` is even a bit faster.