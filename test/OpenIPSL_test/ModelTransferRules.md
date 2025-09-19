# Rules for transforming OpenIPSL models to PowerDynamics.jl

- variable and parameter names should stay consistent
- always look at allready implemented files to get an idea of the style

## Terminal conventions
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

## Parameterstructure and Initialization
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

### Default and Guess values
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

### Parameter Categorization
All **public paramers** (also those lifted from the extended models) go to the @parameters blocks.
All **protect parameters** should end up in a plain `begin...end` block before the equations.

The public parameters a categorized in "fixed" and "free" parameters. Fixed parameters have a default value, free parameters have a guess value.

The protecte parameters fall in different categorations:
  - "derived" parameters: those are used in the equations but are essentialy just shorthands for some simple terms combining multiple other parameters. Those can be placed in plain `begin..end` blocks before the equations (there you essentially define terms `p_derived = p1+p2` will be treated as the term `p1 + p2` in the equation block).
  - "absolute values": sometimes, there are absolute values like heuristic values, those can be treated like "derived" parameters.
  - "initialization parameters": those parameters are only used as intermediate results for initialization. they do not appear in the equations, instead they are only used as `start` properties for variables or in the default equations for "free" parameters. Those should be copied
  to the `begin..end` block as well, but they should be commented out with a not that they are only used to initialize some other variables.
  Intermediate results, which are just used in other commented out init equatiosn but not in the actual dynamic equation should be commented out as well.

## If/Else:
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
