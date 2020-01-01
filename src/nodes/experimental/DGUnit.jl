@doc doc"""
```Julia
DGUnit(;I_r)
```


"""
@DynamicNode DGUnit(Pref, Qref) begin
    MassMatrix()
end  begin
end [[Pmeas, dPmeas], [Qmeas, dQmeas]] begin
    # saturation
    # time delay

    Va_amp = abs(u)
    Va_phase = angle(u)
    I_r = 2//3 * conj(Pref + im * Qref) / Va_amp
    du = i - I_r
    dPmeas = 0.
end

export DGUnit
