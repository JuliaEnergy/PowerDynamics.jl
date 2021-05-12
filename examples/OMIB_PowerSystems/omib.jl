using PowerDynamics
using PowerSystems

raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
system = System(raw_file, dyr_file, runchecks = true)

gen = collect(get_components(DynamicGenerator, system))[1]

ionode = MetaGenerator(gen, verbose=false);

ionode.block
#=
IOBlock :##MGenerator#501 with 4 eqs
  ├ inputs:  i_r(t), i_i(t)
  ├ outputs: u_r(t), u_i(t)
  ├ istates: δ(t), ω(t)
  └ iparams: shaft₊Ω_b, …7…,machine₊e_q
=#

ionode.parameter_names .=> ionode.parameters
#=
9-element Vector{Pair{Symbol, Float64}}:
 :machine₊X_d => 0.2995
   :machine₊R => 0.0
   :shaft₊Ω_b => 1.0
     :mover₊η => 1.0
     :shaft₊D => 2.0
 :mover₊P_ref => 1.0
     :shaft₊H => 3.148
 :machine₊e_q => 1.0
 :shaft₊ω_ref => 50.0
=#

ionode.block.system.eqs
#=
4-element Vector{Equation}:
 Differential(t)(δ(t)) ~ shaft₊Ω_b*(ω(t) - shaft₊ω_ref)
 Differential(t)(ω(t)) ~ (1//2)*(mover₊P_ref*mover₊η - ((i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t)))*(machine₊R*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t))) - (machine₊X_d*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t))))) - (machine₊R*machine₊e_q*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t))))) + machine₊X_d*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t)))*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t))))) - (shaft₊D*(ω(t) - shaft₊ω_ref)))*(shaft₊H^-1)
 u_r(t) ~ sin(-δ(t))*(machine₊X_d*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t))) - (machine₊R*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t)))))) - (cos(-δ(t))*(-machine₊X_d*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t)))) - (machine₊R*machine₊e_q*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t))))))
 u_i(t) ~ cos(-δ(t))*(machine₊X_d*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t))) - (machine₊R*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t)))))) + sin(-δ(t))*(-machine₊X_d*(i_r(t)*sin(δ(t)) - (i_i(t)*cos(δ(t)))) - (machine₊R*machine₊e_q*(i_i(t)*sin(δ(t)) + i_r(t)*cos(δ(t)))))
=#
