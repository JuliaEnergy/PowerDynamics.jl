
@DynamicNode DirectDrivePTO(c, R, L) begin
    @assert c > 0 "Damping coefficient (c) should be >0"
    @assert R >= 0 "Resistance (R) should be >=0"
    @assert L > 0 "Inductance (L) should be >0"
    MassMatrix(m_int=[true]) # Indicates that 'i' has an associated mass
end [[i, di]] begin
    velocity_df = PowerDynamics.velocity_df

    # Ensure the velocity DataFrame is not empty
    @assert !isempty(velocity_df) "Velocity DataFrame should not be empty"

    # Fetch the closest time index and corresponding velocity value
    current_time = PowerDynamics.ts
    closest_time_index = argmin(abs.(velocity_df.Time .- current_time))
    velocity = velocity_df[closest_time_index, :Velocity]
    
end

export DirectDrivePTO


#https://wec-sim.github.io/WEC-Sim/dev/user/advanced_features.html#pto-pto-sim

# RM3 with direct drive linear generator