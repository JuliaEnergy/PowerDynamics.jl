using PowerDynamics
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting

using PowerDynamics.Library
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

using CSV
using DataFrames
using CairoMakie

ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENSAL","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)
names(ref)

# bus 1 is provided from outside
# GENSAL generator parameters from OpenIPSL test
GENSAL_BUS = let
    # GENSAL generator parameters from OpenIPSL test case
    S_b = 100e6
    H = 4.28
    M_b = 100e6
    D = 0
    # Machine reactances
    Xd = 1.84
    Xq = 1.75
    Xpd = 0.41
    Xppd = 0.2
    Xppq = 0.2
    Xl = 0.12
    # Time constants
    Tpd0 = 5.0
    Tppd0 = 0.07
    Tppq0 = 0.09
    # Saturation
    S10 = 0.11
    S12 = 0.39
    # Armature resistance
    R_a = 0.0
    # V_b = 400e3
    # ω_b = 2π*50

    # powerflow results, used to set up pfmodel
    # P_0 = 40e6
    # Q_0 = 5.416582e6
    v_0 = 1.0
    angle_0 = 0.070492225331847

    @named gensal = PSSE_GENSAL(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xppd, Xppq, Xl,
        Tpd0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=false,
        efd_input=false,
    )
    busmodel = MTKBus(gensal; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end


sol = OpenIPSL_SMIB(GENSAL_BUS);

## perform tests for specified variables
# Machine currents
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊iq), "gENSAL.iq") < 7e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊id), "gENSAL.id") < 1.3e-3

# State variables (3 states for salient pole)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊PSIkd), "gENSAL.PSIkd") < 6e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊Epq), "gENSAL.Epq") < 2.5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊PSIppq), "gENSAL.PSIppq") < 1e-3

# Rotor dynamics
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊delta), "gENSAL.delta") < 1.8e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊w), "gENSAL.w") < 8e-6

# Observables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊P), "gENSAL.P") < 1.3e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊Q), "gENSAL.Q") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊anglev), "gENSAL.anglev") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensal₊Vt), "gENSAL.Vt") < 5e-5
