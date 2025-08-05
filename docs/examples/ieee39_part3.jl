#=
# [IEEE39 Bus Tutorial - Part III: Dynamic Simulation](@id ieee39-part3)
=#

include("ieee39_part1.jl")
## add additional init formulas as described in Part II
formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)
## find the initial state
s0 = initialize_from_pf!(nw; verbose=false)
nothing #hide
