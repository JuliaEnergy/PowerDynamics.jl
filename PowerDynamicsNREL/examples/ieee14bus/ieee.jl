using PowerDynamicsNREL
import PowerSimulationsDynamics as PSID
using PowerSystems
import Sundials
import Plots

raw_file = joinpath(@__DIR__, "14bus.raw")
dyr_file = joinpath(@__DIR__, "dyn_data.dyr")
sys = System(raw_file, dyr_file, runchecks = true)

time_span = (0.0, 30.0)
perturbation_trip = PSID.BranchTrip(1.0, "BUS 02-BUS 04-i_4")
sim = PSID.Simulation(PSID.ImplicitModel, sys, pwd(), time_span, perturbation_trip)

####
#### Simulate using PowerSimulationDynamics
####
@time PSID.execute!(
    sim, #simulation structure
    Sundials.IDA(), #Sundials DAE Solver
    dtmax = 0.02,
    reset_simulation = true
); #Arguments: Maximum timestep allowed

plt = Plots.plot(xlabel="time", ylabel="Voltage angle [rad]")
for i in 1:14
    θ = PSID.get_voltage_angle_series(sim, i)
    Plots.plot!(plt, θ, label="θ at Bus $i")
end
plt

####
#### Simulate using PowerDynamics
####
# we piggyback on the initialisation routine of powersystems here
# first we generate an equivalent PowerGrid from the system
# we load `sim.sys` instead of `omib_sys` because it contains initialized machine parameters
pg = PowerGrid(sim.sys);

# we can get the initial condition / operation point from the initialized sim object as well
ic = PSID.get_initial_conditions(sim)
x0 = convert_initial_conditions(pg, sim.sys, ic)

dx = zeros(length(x0))
rhs(pg)(dx, x0, nothing, 0.0)
dx = nd_syms = syms_containing(rhs(pg), "") .=> dx
sort!(dx, by=x->abs(last(x)), rev=true)

fault = PartialLineFailure("BUS 02 -> BUS 04", 0.0, (1.0, time_span[2]))
sol = simulate(fault, pg, x0, time_span; dtmax=0.01);
# sol = simulate(fault, pg, x0, time_span; dtmax=0.01, initializealg=ShampineCollocationInit());

#=
TODO: look at transformer rates? and lines defined in both directions
prob experiment arount with small omib with trafo rates and bidirectional lines
=#
