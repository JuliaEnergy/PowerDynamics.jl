using PowerDynamicsNREL
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

machine = PowerDynamicsNREL.base_machine()

machine_p = [machine.R => 1.0,
             machine.X_d => 1.0,
             machine.e_q => 1.0]

#=
Simple AVR

Takes a bunch of inputs and creats the field voltage. Whatever this means....
Since the machine does not use the v_f input this block has no effect!

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Simple-AVR-[AVRSimple]
=#
AVR = PowerDynamicsNREL.avr_simple()

AVR_p = [AVR.K => 2.0,
         AVR.v_ref => 3.14]

#=
Fixed PSS

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/pss/#Fixed-PSS-[PSSFixed]
=#
PSS = PowerDynamicsNREL.pss_fixed()

PSS_p = [PSS.v_fix => 1.0]

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

shaft = PowerDynamicsNREL.single_mass_shaft()

shaft_p = [shaft.Ω_b => 1.2,
           shaft.ω_ref => 50,
           shaft.D => 4,
           shaft.H => 2.2]

#=
prime mover
we'll just use the fixed version
=#
mover = PowerDynamicsNREL.tg_fixed()

mover_p = [mover.P_ref => 404,
           mover.η => 1.0]

#=
Let's plug it together, shall we??
=#
# the AVR eqs won't show up in the final system therefor we must not
# provide thos in the parameter dictionary
para = Dict(vcat(mover_p, shaft_p, machine_p, AVR_p, PSS_p))
block, p = MetaGenerator(para, mover, shaft, machine, AVR, PSS, verbose=false);
