# simulate OMIB example with PowerSimulatonDynamics

using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
gr()

raw_file = joinpath(@__DIR__, "OMIB.raw")
dyr_file = joinpath(@__DIR__, "OMIB.dyr")
omib_sys = System(raw_file, dyr_file, runchecks = true)

time_span = (0.0, 30.0)

perturbation_trip = BranchTrip(1.0, "BUS 1-BUS 2-i_1")
sim = PSID.Simulation(PSID.ImplicitModel, omib_sys, pwd(), time_span, perturbation_trip)

#V_source_change = SourceBusVoltageChange(1.0, case_source, PSID.θ_source_index, 0.1)
#sim = PSID.Simulation(PSID.ImplicitModel, omib_sys, pwd(), time_span, V_source_change)

print_device_states(sim)
x0_init = PSID.get_initial_conditions(sim)

PSID.execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02,
    reset_simulation = true
); #Arguments: Maximum timestep allowed

angle = get_state_series(sim, ("generator-102-1", :δ));
Plots.plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "rotor angle")

volt = get_voltage_magnitude_series(sim, 102);
Plots.plot(volt, xlabel = "time", ylabel = "Voltage [pu]", label = "V_2")

vangle = get_voltage_angle_series(sim, 102);
Plots.plot(vangle, xlabel = "time", ylabel = "Voltage angle [rad]", label = "Va_2")