using PowerSystems
import PowerSimulationsDynamics as PSID
import Sundials
using PowerDynamicsNREL
using LightGraphs
import Plots

raw_file = joinpath(@__DIR__, "14bus.raw")
dyr_file = joinpath(@__DIR__, "dyn_data.dyr")
sys = System(raw_file, dyr_file, runchecks = true)

time_span = (0.0, 30.0)
perturbation_trip = PSID.BranchTrip(1.0, "BUS 02-BUS 04-i_4")
@time sim = PSID.Simulation(PSID.ImplicitModel, sys, pwd(), time_span, perturbation_trip)

x0_init = PSID.get_initial_conditions(sim)

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

powergrid = PowerGrid(sys; verbose=true)

op = find_operationpoint(powergrid);
timespan = (0.0, 3.0)

fault = PartialLineFailure("BUS 01 -> BUS 02", 1.0, (1.0,2.0))
sol = simulate(fault, powergrid, op.vec, timespan);


#=
Topological elements
- Arc
- Area

Devices
- Branch
  - ACBranch
    ├─ DynamicBranch
    ├─ Line
    ├─ MonitoredLine
    ├─ PhaseShiftingTransformer
    ├─ TapTransformer
    └─ Transformer2W
  - DCBranch
    ├─ HVDCLine
    └─ VSCDCLine
- StaticInjection
    ├─ ElectricLoad
    │  ├─ FixedAdmittance
    │  ├─ ControllableLoad
    │  │  └─ InterruptibleLoad
    │  └─ StaticLoad
    │     └─ PowerLoad
    ├─ Generator
    │  ├─ HydroGen
    │  │  ├─ HydroDispatch
    │  │  ├─ HydroEnergyReservoir
    │  │  └─ HydroPumpedStorage
    │  ├─ RenewableGen
    │  │  ├─ RenewableDispatch
    │  │  └─ RenewableFix
    │  └─ ThermalGen
    │     ├─ ThermalMultiStart
    │     └─ ThermalStandard
    ├─ Source
    ├─ StaticInjectionSubsystem
    │  └─ HybridSystem
    └─ Storage
       ├─ BatteryEMS
       └─ GenericBattery
- DynamicInjection
  ├─ DynamicGenerator
  ├─ DynamicInverter
  └─ PeriodicVariableSource
- RegulationDevice

AbstractTrees.print_tree(Area)
AbstractTrees.print_tree(AggregationTopology)

AbstractTrees.print_tree(Topology)

supertype(LoadZone)


AbstractTrees.print_tree(Bus)
supertype(Bus)

=#
