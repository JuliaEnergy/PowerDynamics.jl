#=
# Modular Inverter Framework
This section describes the modular inverter framework developed in the MARiE-Project funded by the
*German Federal Ministry for Economic Affairs and Climate Action*
=#

using PowerDynamics
using PowerDynamics: IOComponents
using PowerDynamics.ModularInverter


IOComponents.Power

MIParameters

# Foobar

VoltageControlFOM()
VoltageControlCSFOM()

CurrentControlFOM()

DroopControl()
Synchronverter()
FixedVoltage()

PLLCurrent()
ConstantPower()

FiltConstantPower()


balresid(VoltageControlFOM(), 20)
