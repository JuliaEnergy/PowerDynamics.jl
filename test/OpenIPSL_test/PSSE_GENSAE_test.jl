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
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENSAE","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
# GENSAE generator parameters from OpenIPSL test
GENSAE_BUS = let
    # GENSAE generator parameters from OpenIPSL test case
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

    @named gensae = PSSE_GENSAE(;
        S_b, H, M_b, D,
        Xd, Xq, Xpd, Xppd, Xppq, Xl,
        Tpd0, Tppd0, Tppq0,
        S10, S12, R_a,
        pmech_input=false,
        efd_input=false,
    )
    busmodel = MTKBus(gensae; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_SMIB(GENSAE_BUS);

## perform tests for specified variables
# Machine currents
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊iq), "gENSAE.iq") < 7e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊id), "gENSAE.id") < 1.3e-3

# State variables (3 states for salient pole)
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊PSIkd), "gENSAE.PSIkd") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Epq), "gENSAE.Epq") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊PSIppq), "gENSAE.PSIppq") < 1e-3

# Rotor dynamics
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊delta), "gENSAE.delta") < 1.7e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊w), "gENSAE.w") < 8e-6

# Observables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊P), "gENSAE.P") < 1.3e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Q), "gENSAE.Q") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊anglev), "gENSAE.anglev") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :gensae₊Vt), "gENSAE.Vt") < 4e-5
