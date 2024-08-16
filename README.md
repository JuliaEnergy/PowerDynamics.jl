# PowerDynamicsPrototype

[![Build Status](https://github.com/hexaeder/PowerDynamicsPrototype.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hexaeder/PowerDynamicsPrototype.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Basic model

I think the core functionality of PowerDynamics should be to provide a "Bus" type/constructor.

Each bus, contains a list of devices, devices can be probably pure MTK models. They should be either qualified as "current sources" or "voltage sources".

The bus itself can posess meta fields like the "static model" (PV, PQ, ...)

In order to generate the dynamic models there are essentialy 3 ways:
- error if more then 2 devices are voltage sources
- if one device is voltage source
  - voltage source device 
  - current source devices are added like power lines
- if all devices are current sources
  - bus voltage is implicitly give as algebraic constraint
  - alternatively: EMT bus: bus voltage is generated over capacity
