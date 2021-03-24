using PowerDynamics
using ModelingToolkit
using BlockSystems

#=
Base Machine

Inputs: current
Opututs: voltage + electrical torque

i_p,i_q  u_p,u_q
   ↓        ↑
 +---------------+
 | r_a, x_d, e_q |-→ τ_e
 +---------------+

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/machines/#Classical-Model-(Zero-Order)-[BaseMachine]
=#
@parameters t
@parameters r_a x_d e_q i_d(t) i_q(t)
@variables v_d(t) v_q(t) τ_e(t)

machine = IOBlock([v_d ~ -r_a*i_d + x_d*i_q,
                   v_q ~ -x_d*i_d - r_a*i_q * e_q,
                   τ_e ~ (v_q + r_a*i_q)*i_q + (v_d + r_a*i_d)*i_d],
                  [i_d, i_q], [v_d, v_q, τ_e], name=:machine)

machine_p = [machine.r_a => 1.0,
             machine.x_d => 1.0,
             machine.e_q => 1.0]

#=
Simple AVR

Takes a bunch of inputs and creats the field voltage. Whatever this means....
Since the machine does not use the v_f input this block has no effect!

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Simple-AVR-[AVRSimple]
=#
@parameters K v_h(t) v_ref
@variables v_f(t)

controller = IOBlock([v_f ~ K*(v_ref - v_h)],
                     [v_h], [v_f], name=:AVR)

controller_p = [controller.K => 2.0,
                controller.v_ref => 3.14]

#=
Single mass shaft

The mechanical and electric torque acts on the shaft. The shaft
has angular velocity and acceleration.

      +------------------+
τ_e -→| Ω_b, ω_ref, D, H |-→ δ
τ_m -→|                  |-→ ω
      +------------------+

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/shafts/#Rotor-Mass-Shaft-[SingleMass]
=#

dt = Differential(t)
@parameters τ_m(t) τ_e(t) D H ω_ref Ω_b
@variables δ(t) ω(t)

shaft = IOBlock([dt(δ) ~ Ω_b * (ω - ω_ref),
                 dt(ω) ~ 1/(2H) * (τ_m - τ_e - D*(ω - ω_ref))],
                [τ_m, τ_e], [δ, ω], name=:shaft)

shaft_p = [shaft.Ω_b => 1.2,
           shaft.ω_ref => 50,
           shaft.D => 4,
           shaft.H => 2.2]

#=
prime mover
we'll just use the fixed version
=#
@parameters P_ref
@variables τ_m(t)

mover = IOBlock([τ_m ~ P_ref], [], [τ_m], name=:mover)

mover_p = [mover.P_ref => 404]

#=
Let's plug it together, shall we??
=#
# the controller eqs won't show up in the final system therefor we must not
# provide thos in the parameter dictionary
para = Dict(vcat(mover_p, shaft_p, machine_p))

node = MetaGenerator(mover, shaft, machine, controller, para);

symbolsof(node) # u_r, u_i, shaft₊δ shaft₊ω

# let's try this node with totally legit parameters!
using OrderedCollections: OrderedDict
buses=OrderedDict(
    "bus1"=> node,
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches);
operationpoint = find_operationpoint(powergrid);
timespan= (0.0,0.1)

# simulating a voltage perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)
