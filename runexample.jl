using PowerDynamics: read_powergrid, Json, find_operationpoint, rhs, simulate, PowerDrop

powergrid = read_powergrid("ieee-14-minimal.json", Json)

operationpoint = find_operationpoint(powergrid)

result = simulate(PowerDrop(
    fraction = 0.9,
    node_number = 1,
    tspan_fault = (2.,3.)),
    powergrid, operationpoint, (0., 5.))
