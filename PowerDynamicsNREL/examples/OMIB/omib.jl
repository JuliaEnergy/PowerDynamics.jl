using PowerDynamicsNREL
using PowerSystems
using ModelingToolkit
using LinearAlgebra: norm
using Plots
using OrderedCollections: OrderedDict

"take real and imaginary part and return norm and argument"
to_exp(re, im) = norm([re, im]), atan(im, re)
"tage norm an argument and return real and imag part"
to_cart(mag, arg) = mag*cos(arg), mag*sin(arg)

# system = build_system(PSSETestSystems, "psse_OMIB_sys")
raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
system = System(raw_file, dyr_file, runchecks = true)

####
#### How to access the variables
####
collect(get_components(Generator, system))[1]
collect(get_components(DynamicGenerator, system))[1].machine.Xd_p
collect(get_components(Arc, system))[1]
collect(get_components(Branch, system))[1]
collect(get_components(Branch, system))[2]
collect(get_components(ElectricLoad, system))
collect(get_components(Storage, system))
collect(get_components(LoadZone, system))[1]
collect(get_components(Bus, system))[1]
collect(get_components(Bus, system))[2]


# get the dynamic generator
# FIXME: dynamic parameters should be initialized first using the power flow solution
dyngen = deepcopy(collect(get_components(DynamicGenerator, system))[1])

# create a MetaGenerator block and use it to create the BusNode
metagen = MetaGenerator(dyngen; verbose=false);
metagen[2][metagen[1].machine₊e_q] = 1.0278632050189886 # taken from actual initialization with PSD
busnode = BusNode(metagen; name=:bus1);

# symbolsof(busnode)
# busnode.block
# busnode.parameter_names .=> busnode.parameters
# busnode.block.system.eqs

####
#### Simulate the system
####
buses=OrderedDict(
    "bus1"=> SlackAlgebraic(U=1.05),
    "bus2"=> busnode)

branches=OrderedDict(
    "branch1"=> PiModelLine(from="bus1", to="bus2",y=2/(0.1im), y_shunt_km=0, y_shunt_mk=0))

powergrid = PowerGrid(buses, branches);

op = find_operationpoint(powergrid);

####
#### Check the initial state (for comparison with PowerSimulatonDynamics)
####
to_exp(op.vec[1:2]...) # v slack in exp
to_exp(op.vec[3:4]...) # v gen in exp
op.vec[5:6] # δ and ω of shaft

# actual simulation
timespan = (0.0, 100.0)
fault = PartialLineFailure("branch1", 0.5, (1.0,100.0))
sol = simulate(fault, powergrid, op.vec, timespan; dtmax=0.01);

####
#### Check the final state (for comparison with PowerSimulatonDynamics)
####
finalstate = sol.dqsol(99)
to_exp(finalstate[1:2]...) # v slack in exp
to_exp(finalstate[3:4]...) # v gen in exp
finalstate[5:6] # δ and ω of shaft

####
#### Plot the solution
####
plot(sol.dqsol; vars=[:generator_102_1₊δ_2], ylimits=(0.325,0.45))
plot(sol.dqsol; vars=[:generator_102_1₊ω_2])

u1exp = to_exp.(sol.dqsol[:u_r_1],  sol.dqsol[:u_i_1])
u2exp = to_exp.(sol.dqsol[:u_r_2],  sol.dqsol[:u_i_2])
plot(sol.dqsol.t, getindex.(u1exp, 1), label="slack u magnitude")
plot!(sol.dqsol.t, getindex.(u2exp, 1), label="gen u magnitude")
plot(sol.dqsol.t, getindex.(u1exp, 2), label="slack u angle")
plot!(sol.dqsol.t, getindex.(u2exp, 2), label="gen u angle")
