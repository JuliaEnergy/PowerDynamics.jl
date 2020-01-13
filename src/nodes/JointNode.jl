@doc doc"""
```Julia
JointNode(;current_sources, voltage_source)
```

A node type that keeps the current fixed as a desired values `I_r`.

`CSIMinimal` models an inverters as an ideal current source. This can be the
most simple representation of an inverter in grid-feeding mode, according to
Rocabert, Joan, et al. "Control of power converters in AC microgrids." (2012).
Here, additionally to `u`, there are no internal dynamic variables.


# Keyword Arguments
- `I_r`: reference/ desired current


# Mathematical Representation
Using `CSIMinimal` for node ``a`` gives:
```math
0 = I_{r,a} - \left\|i_a\right\|
```

"""
@DynamicNode CSIMinimal(I_r) begin
    MassMatrix()
end  begin
end [] begin
    du = i - I_r
end

export CSIMinimal
