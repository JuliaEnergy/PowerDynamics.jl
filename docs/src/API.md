# API

The following functions are designed for public use.

## Modeling Tools

### Connectors and Base Components
```@docs
Terminal
BusBar
LineEnd
```

### Bus and Line Construction
```@docs
MTKBus
MTKLine
CompositeInjector
```

### Base Unit Calculations
```@docs
Ibase
Zbase
Ybase
```

## Network Components
```@docs
compile_bus
compile_line
simplify_mtkbus
simplify_mtkline
Bus
Line
```

## Interface Checking Functions
```@docs
isinjectormodel
isbusmodel
isbranchmodel
islinemodel
```

## Power Flow Analysis

### Power Flow Bus Types
```@docs
pfSlack
pfPV
pfPQ
```

### Power Flow Solution Functions
```@docs
solve_powerflow
initialize_from_pf!
initialize_from_pf
show_powerflow
powerflow_model
ispfmodel
```

### Power Flow Model Management Functions
```@docs
has_pfmodel
get_pfmodel
set_pfmodel!
delete_pfmodel!
```

## Power Flow Initialization Constraints

### Constraint Types
```@docs
PFInitConstraint
@pfinitconstraint
PFInitFormula
@pfinitformula
```

### Constraint Management Functions
```@docs
add_pfinitconstraint!
add_pfinitformula!
set_pfinitconstraint!
set_pfinitformula!
has_pfinitconstraint
has_pfinitformula
get_pfinitconstraints
get_pfinitformulas
delete_pfinitconstraints!
delete_pfinitformulas!
copy_pf_parameters
```

## Utils
```@docs
refine_timeseries
CallbackVerbosity
set_callback_verbosity!
get_callback_verbosity
with_callback_verbosity
SaturationConfig
SaturationConfiguration
set_saturation_config!
get_saturation_config
with_saturation_config
```
