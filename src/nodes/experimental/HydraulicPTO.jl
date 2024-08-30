
@DynamicNode HydraulicPTO(βe, Vo, Ap, Cd, D, Jt, bg, bf, R, L, Av, rho, alpha, k1) begin
    @assert βe > 0 "Effective bulk modulus (βe) should be >0"
    @assert Vo > 0 "Initial volume (Vo) should be >0"
    @assert Ap > 0 "Piston area (Ap) should be >0"
    @assert Cd > 0 "Discharge coefficient (Cd) should be >0"
    @assert D > 0 "Motor displacement (D) should be >0"
    @assert Jt > 0 "Mass moment of inertia (Jt) should be >0"
    @assert bg >= 0 "Generator damping (bg) should be >=0"
    @assert bf >= 0 "Frictional damping (bf) should be >=0"
    @assert Av > 0 "Orifice area (Av) should be >0"
    @assert rho > 0 "Fluid density (rho) should be >0"
    @assert alpha >= 0 "Swashplate angle ratio (alpha) should be >=0"
    MassMatrix(m_int=[true]) # Indicates that 'i' has an associated mass
end [[i, di, omega]] begin
    velocity_df = PowerDynamics.velocity_df

    # Ensure the velocity DataFrame is not empty
    @assert !isempty(velocity_df) "Velocity DataFrame should not be empty"

    # Fetch the closest time index and corresponding velocity value
    current_time = PowerDynamics.ts
    closest_time_index = argmin(abs.(velocity_df.Time .- current_time))
    velocity = velocity_df[closest_time_index, :Velocity]

    # Access and update global pressure variables
    global initial_pressure_A = PowerDynamics.initial_pressure_A
    global initial_pressure_B = PowerDynamics.initial_pressure_B

    # Pressure rate of change calculations
    pA_dot = (βe / (Vo - Ap * velocity)) * (Ap * velocity)
    pB_dot = (βe / (Vo + Ap * velocity)) * (-Ap * velocity)

    # Assume a fixed time step for simulation
    dt = 0.001  # Time step size; fixed for this model

    # Update global pressures based on pressure rate of change
    PowerDynamics.initial_pressure_A += pA_dot * dt
    PowerDynamics.initial_pressure_B += pB_dot * dt

    # Valve flow calculations
    V1_dot = Cd * Av * sqrt(2 / rho * abs(initial_pressure_A - initial_pressure_B)) * tanh(k1 * (initial_pressure_A - initial_pressure_B))
    V2_dot = Cd * Av * sqrt(2 / rho * abs(initial_pressure_B - initial_pressure_A)) * tanh(k1 * (initial_pressure_B - initial_pressure_A))

    # Accumulator flow calculations
    VC_dot = -alpha * D * PowerDynamics.initial_omega + V1_dot + V2_dot
    VD_dot = alpha * D * PowerDynamics.initial_omega - V1_dot - V2_dot

    # Motor and generator dynamics
    omega_dot = (1 / Jt) * (D * (initial_pressure_A - initial_pressure_B) - bg * PowerDynamics.initial_omega - bf * PowerDynamics.initial_omega)

    # Update global omega based on omega_dot
    PowerDynamics.initial_omega += omega_dot * dt

    # PTO force
    F_pto = (initial_pressure_A - initial_pressure_B) * Ap

    # Electromotive force (emf) and di/dt
    emf = F_pto
    di = (emf - R * i[1]) / L
    
    # Voltage (du) across the PTO
    du = emf - R * i[1] - L * di
end

export HydraulicPTO
