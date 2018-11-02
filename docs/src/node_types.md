# Node Types

The currently implementes node types are

* Purely Algebraic:
  * [`PowerDynBase.PQAlgebraic`](@ref) (PQ-bus)
  * [`PowerDynBase.PVAlgebraic`](@ref) (PV-bus)
  * [`PowerDynBase.SlackAlgebraic`](@ref) (Slack-bus / Vφ-bus)
* Synchronous Machine Models:
  * [`PowerDynBase.SwingEq`](@ref) (2nd order)
  * [`PowerDynBase.SwingEqLVS`](@ref) (2nd order with an additional term for numerical voltage stability)
  * [`PowerDynBase.FourthEq`](@ref) (4th order)
* Voltage Source Inverters:
  * [`PowerDynBase.VSIMinimal`](@ref)
  * [`PowerDynBase.VSIVoltagePT1`](@ref)


```@docs
PQAlgebraic
PVAlgebraic
SlackAlgebraic
SwingEq
SwingEqLVS
FourthEq
VSIMinimal
VSIVoltagePT1
```
