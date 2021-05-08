using PowerDynamics
using PowerSystems

raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
system = System(raw_file, dyr_file, runchecks = true)

gen = collect(get_components(DynamicGenerator, system))[1]

ionode = MetaGenerator(gen, verbose=false);

ionode.block
#=
IOBlock :##MGenerator#257 with 4 eqs
  ├ inputs:  i_r(t), i_i(t)
  ├ outputs: u_r(t), u_i(t)
  ├ istates: shaft₊δ(t), shaft₊ω(t)
  └ iparams: shaft₊ω_ref, …7…,machine₊e_q
=#

ionode.parameter_names .=> ionode.parameters
#=
9-element Vector{Pair{Symbol, Float64}}:
 :shaft₊ω_ref => 50.0
 :machine₊X_d => 0.2995
   :machine₊R => 0.0
 :machine₊e_q => 1.0
     :mover₊η => 1.0
   :shaft₊Ω_b => 1.0
 :mover₊P_ref => 1.0
     :shaft₊D => 2.0
     :shaft₊H => 3.148
=#

ionode.block.system.eqs
#=
Differential(t)(shaft₊δ(t)) ~ shaft₊Ω_b*(shaft₊ω(t) - shaft₊ω_ref)
Differential(t)(shaft₊ω(t)) ~ (1//2)*(mover₊P_ref*mover₊η - ((i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t)))*(machine₊R*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t))) - (machine₊X_d*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t))))) - (machine₊R*machine₊e_q*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t))))) + machine₊X_d*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t)))*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t))))) - (shaft₊D*(shaft₊ω(t) - shaft₊ω_ref)))*(shaft₊H^-1)
u_r(t) ~ sin(-shaft₊δ(t))*(machine₊X_d*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t))) - (machine₊R*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t)))))) - (cos(-shaft₊δ(t))*(-machine₊X_d*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t)))) - (machine₊R*machine₊e_q*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t))))))
u_i(t) ~ cos(-shaft₊δ(t))*(machine₊X_d*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t))) - (machine₊R*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t)))))) + sin(-shaft₊δ(t))*(-machine₊X_d*(i_r(t)*sin(shaft₊δ(t)) - (i_i(t)*cos(shaft₊δ(t)))) - (machine₊R*machine₊e_q*(i_i(t)*sin(shaft₊δ(t)) + i_r(t)*cos(shaft₊δ(t)))))
=#
