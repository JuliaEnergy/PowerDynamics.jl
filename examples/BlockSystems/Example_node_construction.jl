
#=
## Example Usage of MTK components
This example shall illustrate how to build a node from MTK components. For this purpose a gridforming inverter was chosen with the below block diagram.
It contains a low-pass filters and a droop control component for active and reactive power, respectively.
Since the inverter is grid-forming it also contains a voltage source 
Implementation of `VSIVoltagePT1` as an IONode.

```

                +----------+ +---------+
        +-----+ | p_filter |-| p_droop |  +----------+
    i_r-|    Q|-|   τ_P    | | P, ω_r  |--|          |
    i_i-|     | +----------+ +---------+ ω| v_source |--+-u_r
        | pow |                           |   τ_v    |  |
     +--|     | +----------+ +---------+ V|          |-+|-u_i
     |+-|    P|-| q_filter |-| q_droop |--|          | ||
     || +-----+ |   τ_Q    | | Q, V_r  |  +----------+ ||
     ||         +----------+ +---------+               ||
     |+------------------------------------------------+|
     +----

```
=#

using BlockSystems
using PowerDynamics.IOComponents
using PowerDynamics: IONode, GFI, SlackAlgebraic, PiModelLine, PowerGrid, find_operationpoint, ChangeInitialConditions, simulate, Inc
using OrderedCollections: OrderedDict
#=
The parameters for initalization are described in the comments below
GFI(;τ_v,τ_P,τ_Q,K_P,K_Q,V_r,P,Q,ω_r)
=#
name = gensym(:GridForming)
parameters = Dict(
        :τ_v => 0.005,    # time constant voltage control delay
        :τ_P => 0.5,      # time constant active power measurement
        :τ_Q => 0.5,      # time constant reactive power measurement
        :K_P => 0.396,     # droop constant frequency droop
        :K_Q => 0.198,     # droop constant voltage droop
        :V_r => 1.0,   # reference/ desired voltage
        :P   => 0.303, # active (real) power infeed
        :Q   => 0.126, # reactive (imag) power infeed                .
        :ω_r => 0.0)   # refrence/ desired frequency

#=
In the following components that build the grid-forming inverter are defined with implementions from the IOComponents library. They are defined with a 
name and the necessary parameters for initalization.
=#
p_filter = LowPassFilter(name=:p_filter, τ=:τ_P)
q_filter = LowPassFilter(name=:q_filter, τ=:τ_Q)
q_droop = DroopControl(name=:q_droop, x_ref=:Q, u_ref=:V_r, K=:K_Q)
v_source = VoltageSource(name=:v_source, τ=:τ_v)

#=
Only difference to what is usually part of a block diagram is the power component that calculates active and reactive power 
from the grid current and voltage. 
=#
pow = Power(name=:pow)

#=
The active power droop control could be also defined with the DroopControl component. However, since the strength of the modularity
is that you can also define your own components and plug'n'play with the already existing components, we simply define our own droop control.
So in the following a self-defined droop control "MyDroopControl" for active power is implemented.


```
    u = - K*(x - x_ref) + u_ref

           +-----------------+
    x(t) --| K, x_ref, u_ref |-- u(t)
           +-----------------+

```
    IOBlock :##droop# with 1 eqs
    ├ inputs:  x(t)
    ├ outputs: u(t)
    ├ istates: (empty)
    └ iparams: K, x_ref, u_ref
=#

function MyDroopControl(;name=gensym(:droop), renamings...)
    @parameters t K x_ref u_ref
    @parameters x(t)
    @variables u(t)
    D = Differential(t)

    block = IOBlock([u ~ - K * (x - x_ref) + u_ref], # output is the droop voltage v
                    [x], [u]; name)

    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

p_droop = MyDroopControl(name=:p_droop, x_ref=:P, u_ref=:ω_r, K=:K_P)


gfi = IOSystem([pow.P => p_filter.input,
                pow.Q => q_filter.input,
                p_filter.output => p_droop.x,
                q_filter.output => q_droop.x,
                p_droop.u => v_source.ω,
                q_droop.u => v_source.V,
                v_source.u_r => pow.u_r,
                v_source.u_i => pow.u_i],
               [pow, p_filter, p_droop, q_filter, q_droop, v_source],
               name=name,
               outputs=[v_source.u_i, v_source.u_r],
               namespace_map=[p_droop.P => :P,
                              q_droop.Q => :Q])


connected = connect_system(gfi)
IONode(connected, parameters)

#=
Finally we build a small grid example with new node defined with our IOComponents and a Slack bus connected with a pi-model line
=#
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

#=
Then a simple voltage perturbation at node is simulated.
=#
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)

