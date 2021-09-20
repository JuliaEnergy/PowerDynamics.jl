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

machine_b = PowerDynamicsNREL.base_machine()

machine_p = Dict(machine_b.R => 1.0,
                 machine_b.X_d => 1.0,
                 machine_b.e_q => 1.0)

machine = BlockPara(machine_b, machine_p)

#=
Simple AVR

Takes a bunch of inputs and creats the field voltage. Whatever this means....
Since the machine does not use the v_f input this block has no effect!

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Simple-AVR-[AVRSimple]
=#
AVR_b = PowerDynamicsNREL.avr_simple()

AVR_p = Dict(AVR_b.K => 2.0,
             AVR_b.v_ref => 3.14)

AVR = BlockPara(AVR_b, AVR_p)

#=
Fixed PSS

https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/pss/#Fixed-PSS-[PSSFixed]
=#
PSS_b = PowerDynamicsNREL.pss_fixed()

PSS_p = Dict(PSS_b.v_fix => 1.0)

PSS = BlockPara(PSS_b, PSS_p)

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

shaft_b = PowerDynamicsNREL.single_mass_shaft()

shaft_p = Dict(shaft_b.Ω_b => 1.2,
               shaft_b.ω_ref => 50,
               shaft_b.D => 4,
               shaft_b.H => 2.2)

shaft = BlockPara(shaft_b, shaft_p)

#=
prime mover
we'll just use the fixed version
=#
mover_b = PowerDynamicsNREL.tg_fixed()

mover_p = Dict(mover_b.P_ref => 404,
               mover_b.η => 1.0)

mover = BlockPara(mover_b, mover_p)

#=
Let's plug it together, shall we??
=#
# the AVR eqs won't show up in the final system therefor we must not
blockp = MetaGenerator(mover, shaft, machine, AVR, PSS, verbose=false);
