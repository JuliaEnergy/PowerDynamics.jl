# Instructions for porting OpenIPSL models to PowerDynamics.jl

## Test Setup

To implement and validate a new OpenIPSL model in PowerDynamics.jl, you need to create the following file structure and follow this workflow:

### Required Files Structure
```
test/OpenIPSL_test/
├── MODEL_NAME_test.jl           # Main Julia test file
└── MODEL_NAME/                  # Reference data directory
    ├── ref_data_gen.sh          # Script to generate reference data
    ├── model_simulation.mos     # OpenModelica simulation script
    └── modelica_results.csv.gz  # Compressed reference data (generated)

src/Library/OpenIPSL/
└── Category/                    # e.g., Machines/, Loads/, etc.
    └── MODEL_NAME.jl           # PowerDynamics model implementation
```

### Workflow Steps

1. **Create PowerDynamics Model Implementation**
   - Implement the model in `src/Library/OpenIPSL/Category/MODEL_NAME.jl`
   - Follow the transformation rules detailed below
   - Ensure proper parameter handling and terminal conventions

2. **Set Up Reference Data Generation**
   - Create `MODEL_NAME/ref_data_gen.sh` that:
     - Clones OpenIPSL at a specific version for reproducibility
     - Sets `OPENIPSL_PATH` environment variable
     - Runs the OpenModelica simulation script
     - Compresses and stores the results
   - Create `MODEL_NAME/model_simulation.mos` that:
     - Loads the OpenIPSL library
     - Simulates the corresponding OpenIPSL test case
     - Extracts relevant variables to CSV format using `filterSimulationResults`
     - for debuggin, it is sometimes usefull to alxo export a "extended" version, that one is ignored in gitignore and should not ever be added to git.

    **BOTH scripts can by copied from GENCLS tests and slighlty adapted!**

3. **Create Validation Test**
   - Implement `MODEL_NAME_test.jl` that:
     - Loads reference data from the compressed CSV
     - Sets up the PowerDynamics model with equivalent parameters
     - Runs simulation using `OpenIPSL_SMIB()` or appropriate test scenario
     - Compares key variables using `ref_rms_error()` with appropriate tolerances
     - Uses `@test` assertions to validate RMS errors are below thresholds

4. **Integration with Test Suite**
   - Add the test to `test/runtests.jl` under the "OpenIPSL Model Tests" section:
     ```julia
     @safetestset "MODEL_NAME" begin include(joinpath("OpenIPSL_test", "MODEL_NAME_test.jl")) end
     ```

### Key Testing Functions
- `ref_rms_error(sol, ref, VIndex(), "ref_var_name")`: Computes RMS error between simulation and reference
- `VIndex(bus_name, :component₊variable)`: Indexes into solution for specific variables
- `OpenIPSL_SMIB(bus_model)`: Single Machine Infinite Bus test scenario

The validation should ensure RMS errors are below appropriate thresholds (typically 1e-3 to 1e-5 depending on variable type).

The OpenIPSL folder is available locally! No need to web searches, any questions about OpenIPSL are best answered in the defining Modelica files.

## Rules for transforming OpenIPSL models to PowerDynamics.jl

**BEFORE ACTUALLY STARTING, LOOK AT EXISTING PORTS TO UNDERSTAND HOW TO DO IT**

- variable and parameter names should stay consistent
- always look at allready implemented files to get an idea of the style

### Terminal conventions
Both modeling systems use different terminal conventions. To make the equations comparable, introduce the variables
```julia
@variables begin
    # OpenIPSL-style variables for equation compatibility
    pir(t), [guess=0, description="Real part of terminal current (OpenIPSL convention) [pu]"]
    pii(t), [guess=0, description="Imaginary part of terminal current (OpenIPSL convention) [pu]"]
    pvr(t), [guess=1, description="Real part of terminal voltage (OpenIPSL convention) [pu]"]
    pvi(t), [guess=0, description="Imaginary part of terminal voltage (OpenIPSL convention) [pu]"]
end
```
with the equations
```julia
@equations begin

    # Conversion between PowerDynamics and OpenIPSL conventions
    pir ~ -terminal.i_r
    pii ~ -terminal.i_i
    pvr ~ terminal.u_r
    pvi ~ terminal.u_i
end
```
and then use `pir` instead of `p.ir` and so on.

### Parameterstructure and Initialization
The initialization in OpenIPSL and PowerDynamics is quite different: while
OpenIPSL uses a "explicit" Initialization (i.e. every variable explicitly calculated) while
PowerDynamics just uses some guesses and then numericially solves the initialization problem in 
order to find the correct values.

The main culprit here is the handling of parameters. First, you need to collect and categorized all parameters in from the OpenIPSL models. The main distiction to make here is between **public parameters** and **protected parameters**.
Some parameters hide in the inheritence chain, we normally don't want to copy the full extension
structure. Most notably, there is this `pfCompontent`
```
extends OpenIPSL.Electrical.Essentials.pfComponent(
    final enabledisplayPF=false,
    final enableangle_0=true,
    final enablev_0=true,
    final enableQ_0=true,
    final enableP_0=true,
    final enablefn=false,
    final enableV_b=false,
    final enableS_b=false
);
```
which tells you what additonal parameters to include, for example in this case we add three parameters to our list: `v_0`, `angle_0`, `P_0` and `Q_0`. Those parameters are considered
"free" parameters for our initialization. The other parameters don't need to be included.

#### Default and Guess values
In ModelingTolkit/PowerDynamics, variables can have default and guess values:
```
u(t)=<defaultvalue>, [guess=<guessvalue>, description=..., ...]
```
A default value is fixed in initialization, a guess is just a start value for the initialization
solver. Variables don't get default values, but they all get guess values. For parameters it depends, parameteres which are "infered" from the powerflow, such as the ones from pfComponent
don't get default values, they will be figured out during initialization. They should get
guess values however.
Normal "parameters" should copy the default values from OpenIPSL. If the default value from
OpenIPSL is a term rather than an actual value, thos need to be considered free and only get

#### Parameter Categorization
All **public paramers** (also those lifted from the extended models) go to the @parameters blocks.
All **protect parameters** should end up in a plain `begin...end` block before the equations.
**Protected parameters**, which do not show up in equations but only in "start" values of others keep be ignored! We don't need to keep them around as the initialization works differently in PowerDynamics.

The public parameters a categorized in "fixed" and "free" parameters. Fixed parameters have a default value, free parameters have a guess value.

The protecte parameters fall in different categorations:
  - "derived" parameters: those are used in the equations but are essentialy just shorthands for some simple terms combining multiple other parameters. Those can be placed in plain `begin..end` blocks before the equations (there you essentially define terms `p_derived = p1+p2` will be treated as the term `p1 + p2` in the equation block).
  - "absolute values": sometimes, there are absolute values like heuristic values, those can be treated like "derived" parameters.
  - "initialization parameters": those parameters are only used as intermediate results for initialization. they do not appear in the equations, instead they are only used as `start` properties for variables or in the default equations for "free" parameters. Those can be ignored and don't need porting! 

### If/Else:
In ModelingToolkit, you cannot use if / else statmenes in eqautions. Instead, we need the fucntional form:

```julia
## WRONG!!
@equations begin
    if a < 0 && b > 1
        vf ~ x
    else
        vf ~ x^2
    end
end
```
```julia
## CORRECT!!
@equations begin
    vf ~ ifelse(a<0 & b>1,
        x,
        x^2
    )
end
```

## OpenIPSL Models Available for Implementation

The following OpenIPSL models use the SMIB test harness and are candidates for PowerDynamics.jl implementation. Each model represents a complete test case with reference data generation capability.

### Machine Models (PSSE)
- [X] **GENCLS** - Classical generator model
  - *Dependencies: None (standalone model)*
- [ ] **GENSAL** - Salient pole generator
  - *Dependencies: None (standalone model, connects PMECH0→PMECH and EFD0→EFD)*
- [ ] **GENSAE** - Salient pole generator with saturation
  - *Dependencies: None (standalone model, connects PMECH0→PMECH and EFD0→EFD)*
- [ ] **GENROE** - Round rotor generator with saturation
  - *Dependencies: None (standalone model, connects PMECH0→PMECH and EFD0→EFD)*
- [X] **GENROU** - Round rotor generator
  - *Dependencies: None (standalone model, connects PMECH0→PMECH and EFD0→EFD)*
- [ ] **GEN** - Generic generator model
  - *Dependencies: Uses GENROE as machine, ConstantPower governor, ConstantExcitation exciter, DisabledPSS*

### Excitation Systems (PSSE)
- [ ] **EXST1** - Static excitation system
  - *Dependencies: Requires GENROE machine model*
- [ ] **ST5B** - IEEE Type ST5B exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **ESDC1A** - DC exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **URST5T** - Brushless exciter
  - *Dependencies: Requires GENROU machine model*
- [ ] **SEXS** - Simplified excitation system
  - *Dependencies: Requires GENROU machine model*
- [ ] **SCRX** - Static exciter
  - *Dependencies: Requires GENROU machine model*
- [ ] **ESST1A** - Static exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **ESST4B** - Static exciter
  - *Dependencies: Requires GENROU machine model*
- [ ] **EXAC1** - AC exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **EXAC2** - AC exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **ESDC2A** - DC exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **IEEET2** - IEEE Type 2 exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **ESAC1A** - AC exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **IEEEX1** - IEEE Type 1 exciter
  - *Dependencies: Requires GENROE machine model*
- [ ] **EXNI** - Exciter model
  - *Dependencies: Requires GENROE machine model*
- [ ] **ESAC2A** - AC exciter
  - *Dependencies: Requires GENROU machine model*

### Turbine Governors (PSSE)
- [ ] **GGOV1** - General governor
  - *Dependencies: Requires GENROU machine model*
- [ ] **IEEEG1** - IEEE general governor
  - *Dependencies: Requires GENROU machine model + IEEET1 excitation system*
- [ ] **GAST** - Gas turbine governor
  - *Dependencies: Requires GENROU machine model + IEEET1 excitation system*
- [ ] **IEESGO** - IEEE standard governor
  - *Dependencies: Requires GENSAL machine model + SCRX excitation system*
- [ ] **GGOV1DU** - General governor with droop
  - *Dependencies: Requires GENROU machine model*
- [ ] **TGOV1** - Steam turbine governor
  - *Dependencies: Requires GENROE machine model*
- [ ] **HYGOV** - Hydro turbine governor
  - *Dependencies: Requires GENSAL machine model + SCRX excitation system*

### Power System Stabilizers (PSSE)
- [ ] **PSS2A** - IEEE PSS2A model
  - *Dependencies: Requires GENROE machine model + ESST1A excitation system*
- [ ] **PSS2B** - IEEE PSS2B model
  - *Dependencies: Requires GENROE machine model + ESST1A excitation system*
- [ ] **IEEEST** - IEEE stabilizer
  - *Dependencies: Requires GENROE machine model + ESST1A excitation system*

### CGMES Models
- [ ] **GovHydroIEEE0_Test** - CGMES hydro governor
- [ ] **ExcSEXS** - CGMES simplified exciter

### Solar Models (PSAT)
- [ ] **SolarPVTest** - Solar PV model
- [ ] **SolarPQTest** - Solar PQ model

### Bank Models (PSSE)
- [ ] **CSVGN1** - Static VAR compensator

### Power Plant Examples
- [ ] **Anderson** - Anderson power plant model
- [ ] **IEEE421** - IEEE 421 power plant model


## SMIB Test Harness Derivatives

### SMIBRenewable Test Harness
The SMIBRenewable harness is designed specifically for renewable energy source validation with enhanced measurement capabilities (voltage and current sensors). It uses a GENCLS infinite bus generator.

**Models using SMIBRenewable:**
- [ ] **PVPlant** - PV source with REPCA plant controller model
  - *Dependencies: Requires REGCA1 inverter interface + REECB1 electrical controller + REPCA1 plant controller*
- [ ] **WindPlant** - Wind source with REECCA electrical controller model
  - *Dependencies: Requires REGCA1 inverter interface + REECA1 electrical controller + REPCA1 plant controller + WTDTA1 drive train*
- [ ] **BESSPlant** - Battery Energy Storage System with REPCA plant controller
  - *Dependencies: Requires REGCA1 inverter interface + REECCU1 electrical controller + REPCA1 plant controller*

### SMIBAddOn Test Harness
The SMIBAddOn harness is a specialized version for testing additional renewable energy features with a voltage source instead of a generator for the infinite bus.

**Models using SMIBAddOn:**
- [ ] **PVPlantSolarIrradiance** - PV source with solar irradiance to power conversion
  - *Dependencies: Requires REGCA1 inverter interface + REECB1 electrical controller + REPCA1 plant controller + IrradianceToPower add-on block*
