using PowerDynamics
using ModelingToolkit
using BlockSystems

function getgenerator(;name)
    machine = PowerDynamics.base_machine()
    AVR = PowerDynamics.avr_simple()
    PSS = PowerDynamics.pss_fixed()
    shaft = PowerDynamics.single_mass_shaft()
    mover = PowerDynamics.tg_fixed()

    para = Dict(machine.R => 1.0,
                machine.X_d => 1.0,
                machine.e_q => 1.0,
                AVR.K => 2.0,
                AVR.v_ref => 3.14,
                PSS.v_fix => 1.0,
                shaft.Ω_b => 1.2,
                shaft.ω_ref => 50,
                shaft.D => 4,
                shaft.H => 2.2,
                mover.P_ref => 404,
                mover.η => 1.0)

    MetaGenerator(para, mover, shaft, machine, AVR, PSS; verbose=false, name)
end

@named gen1 = getgenerator()
@named load = BusLoad()
load_p = Dict(load.P => 1.0, load.Q => 0.0)

@named busnode = BusNode(gen1, (load, load_p))

busnode.mass_matrix
busnode.block

symbolsof(busnode)

symbolsof(node) # u_r, u_i, shaft₊δ shaft₊ω

# let's try this node with totally legit parameters!
using OrderedCollections: OrderedDict
buses=OrderedDict(
    "bus1"=> busnode,
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches);
operationpoint = find_operationpoint(powergrid);
timespan= (0.0,0.1)

# simulating a voltage perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
solution1
using Plots
plot(solution1.dqsol)
