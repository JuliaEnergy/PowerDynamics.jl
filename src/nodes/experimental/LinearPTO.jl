
using CSV
using DataFrames
using Interpolations
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
@DynamicNode LinearPTO(
    τ_P, τ_Q, K_P, K_Q, V_r, Q, K_pto, C_pto, η, base_power, scaling_factor, wec_sim_path
) begin
    @assert τ_P > 0 "time constant active power measurement should be >0"
    @assert τ_Q > 0 "time constant reactive power measurement should be >0"
    @assert K_Q > 0 "reactive power droop constant should be >0"
    @assert K_P > 0 "active power droop constant reactive power measurement should be >0"

    # build interpolants objects
    df = CSV.read(wec_sim_path, DataFrame)
    df.Relative_Displacement = df.Float_Position .- df.Spar_Position
    df.Relative_Velocity = df.Float_Velocity .- df.Spar_Velocity
    rv_interp = LinearInterpolation(df.Time, df.Relative_Velocity, extrapolation_bc=Line())
    rd_interp = LinearInterpolation(df.Time, df.Relative_Displacement, extrapolation_bc=Line())

end [[ω, dω]] begin

    relative_velocity = rv_interp(t)
    relative_displacement = rd_interp(t)
    if !isfinite(t)
        relative_velocity = 0.0
        relative_displacement = 0.0
    else
        relative_velocity = rv_interp(t)
        relative_displacement = rd_interp(t)
    end
    # Calculate the power from the PTO 
    F_pto = -K_pto * relative_displacement - C_pto * relative_velocity
    P_mech = -F_pto * relative_velocity
    P_elec = η * P_mech
    P = (P_elec / base_power) * scaling_factor

    p = real(u * conj(i))
    q = imag(u * conj(i))
    dϕ = ω
    v = abs(u)
    dv = 1/τ_Q*(-v + V_r- K_Q *(q-Q))
    du = u * 1im * dϕ + dv*(u/v)
    dω = 1/τ_P*(-ω-K_P*(p-P))

end
