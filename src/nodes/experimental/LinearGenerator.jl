@doc doc"""
```Julia
LinearPermanentMagnetGenerator(;B, l, R, L, velocity)
```

"""
@DynamicNode LinearGenerator(B, l, R, L, velocity) begin
    @assert B > 0 "Magnetic flux density (B) should be >0"
    @assert l > 0 "Effective length (l) should be >0"
    @assert R >= 0 "Resistance (R) should be >=0"
    @assert L > 0 "Inductance (L) should be >0"
    @assert velocity != 0 "Velocity should not be zero"
end [[i, di]] begin
    emf = B * l * velocity
    di = (emf - R * i[1]) / L
    du = emf - R * i[1] - L * di  # enforcing voltage conservation
end

export LinearGenerator
