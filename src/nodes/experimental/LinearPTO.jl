

@doc doc"""
```Julia
LinearPTO(;τ_P,τ_Q,K_P,K_Q,V_r,Q,K_pto, C_pto, η, base_power)
```
A node type that applies frequency and voltage droop control to regulate the dynamics of a wave energy converter's (WEC) linear Power Take-Off (PTO) system.

LinearPTO models the PTO system as an AC voltage source, allowing control over amplitude and frequency (grid-forming mode). The model incorporates time series data from WEC simulations, adjusting power output based on relative velocity and displacement. Power is converted to per-unit using a specified base_power.

this node is based on VSIMinimal

Keyword Arguments
- τ_P: Time constant for active power measurement
- τ_Q: Time constant for reactive power measurement
- K_P: Droop constant for frequency control
- K_Q: Droop constant for voltage control
- V_r: Reference/desired voltage
- Q: Reactive power infeed
- K_pto: Stiffness of the PTO
- C_pto: Damping of the PTO
- η: Efficiency of the PTO system
- base_power: Base power for per-unit conversion

"""
@DynamicNode LinearPTO(τ_P, τ_Q, K_P, K_Q, V_r, Q, K_pto, C_pto, η, base_power, scaling_factor) begin
    @assert τ_P > 0 "time constant active power measurement should be >0"
    @assert τ_Q > 0 "time constant reactive power measurement should be >0"
    @assert K_Q > 0 "reactive power droop constant should be >0"
    @assert K_P > 0 "active power droop constant reactive power measurement should be >0"
end [[ω, dω]] begin

    # handle time series from WEC simulation 
    current_time = PowerDynamics.ts
    closest_time_index = argmin(abs.(wec_simulation_df.Time .- current_time))
    relative_velocity = wec_simulation_df[closest_time_index, :Relative_Velocity]
    relative_displacement = wec_simulation_df[closest_time_index, :Relative_Displacement]

    hs = PowerDynamics.wave_data_df[closest_time_index, :Wave_Height]

    # Calculate the power from the PTO 
    F_pto = -K_pto * relative_displacement - C_pto * relative_velocity
    P_mech = F_pto * relative_velocity
    P_elec = η * P_mech

    P_elec_scaled = scaling_factor * P_elec

    # Convert to per-unit
    P = P_elec / base_power

    p = real(u * conj(i))
    q = imag(u * conj(i))
    dϕ = ω
    v = abs(u)
    dv = 1/τ_Q*(-v + V_r- K_Q *(q-Q))
    du = u * 1im * dϕ + dv*(u/v)
    dω = 1/τ_P*(-ω-K_P*(p-P))

    # @info "P_mech : $P_mech"
    # @info "P_elec : $P_elec"
    # @info "P : $P"
    # @info "ts : $current_time"
    # @info "Time index: $closest_time_index"
    # @info "hs : $hs"
    # @info "Relative Velocity: $relative_velocity"
    # @info "Relative Displacement: $relative_displacement"
    # @info "Call number: $(PowerDynamics.temp)"
    # @info "---- End of Step ----"
    # PowerDynamics.temp += 1
end

export LinearPTO
