# simulate OMIB example with PowerSimulatonDynamics

using PowerDynamicsNREL
import PowerSimulationsDynamics as PSID
using PowerSystems
using Sundials
using OrdinaryDiffEq
using LinearAlgebra
import Plots

raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
omib_sys = System(raw_file, dyr_file, runchecks = true)

time_span = (0.0, 30.0)
perturbation_trip = PSID.BranchTrip(1.0, "BUS 1-BUS 2-i_1")
sim = PSID.Simulation(PSID.ImplicitModel, omib_sys, pwd(), time_span, perturbation_trip)

#V_source_change = SourceBusVoltageChange(1.0, case_source, PSID.θ_source_index, 0.1)
#sim = PSID.Simulation(PSID.ImplicitModel, omib_sys, pwd(), time_span, V_source_change)


####
#### Simulate using PowerSimulationDynamics
####
# PSID.print_device_states(sim)

PSID.execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02,
    reset_simulation = true
); #Arguments: Maximum timestep allowed

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

fault = PartialLineFailure("BUS 1 -> BUS 2", 0.5, (1.0, time_span[2]))
sol = simulate(fault, pg, x0, time_span; dtmax=0.01, initializealg=BrownFullBasicInit());
# sol = simulate(fault, pg, x0, time_span; dtmax=0.01, initializealg=ShampineCollocationInit());

####
#### Plot solutions
####
# rotor angle
angle = PSID.get_state_series(sim, ("generator-102-1", :δ));
Plots.plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "PSD")
Plots.plot!(sol.dqsol; vars=[Symbol("generator-102-1₊δ_2")], label="PD")
Plots.savefig("OMIB_RotorAngle.png")

# voltage plots
"take real and imaginary part and return norm and argument"
to_exp(re, im) = norm([re, im]), atan(im, re)
"tage norm an argument and return real and imag part"
to_cart(mag, arg) = mag*cos(arg), mag*sin(arg)
u1exp = to_exp.(sol.dqsol[:u_r_1],  sol.dqsol[:u_i_1]);
u2exp = to_exp.(sol.dqsol[:u_r_2],  sol.dqsol[:u_i_2]);

# voltage magnitude
slack_mag = PSID.get_voltage_magnitude_series(sim, 101);
gen_mag = PSID.get_voltage_magnitude_series(sim, 102);
Plots.plot(slack_mag, xlabel = "time", ylabel = "Voltage Magnitude [pu]", label="PSD Slack")
Plots.plot!(gen_mag, label="PSD Generator")
Plots.plot!(sol.dqsol.t, getindex.(u1exp, 1), label="PD Slack")
Plots.plot!(sol.dqsol.t, getindex.(u2exp, 1), label="PD Generator")
Plots.savefig("OMIB_VoltageMagnitude.png")

# voltage angle
slack_angle = PSID.get_voltage_angle_series(sim, 101);
gen_angle = PSID.get_voltage_angle_series(sim, 102);
Plots.plot(slack_angle, xlabel = "time", ylabel = "Voltage angle [rad]", label="PSD Slack")
Plots.plot!(gen_angle, label="PSD Generator")
Plots.plot!(sol.dqsol.t, getindex.(u1exp, 2), label="PD Slack")
Plots.plot!(sol.dqsol.t, getindex.(u2exp, 2), label="PD Generator")
Plots.savefig("OMIB_VoltageAngle.png")
