using OpPoDyn
using OpPoDyn.Library
using OpPoDynTesting
using ModelingToolkit
using NetworkDynamics
isinteractive() && using GLMakie

@testset "DynawoPiLine" begin
    @named branchA = DynawoPiLine(XPu=0.022522)
    @named branchB = DynawoPiLine(XPu=0.04189)
    line = Line(LineModel(branchA, branchB));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_1" toi

    @named branchA = DynawoPiLine(XPu=0.022522, RPu=0.01)
    line = Line(LineModel(branchA));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_2" toi
end

@testset "Swing bus" begin
    @named swing = Swing(Pm=1, D=0.1, M=0.005)
    bus = Bus(BusModel(swing));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "SwingBus_1" toi

    # swing bus with load
    @named swing = Swing(Pm=1, D=0.1, M=0.005)
    @named pqload = PQLoad(Pset=-0.5, Qset=-0.2)
    bm = BusModel(swing, pqload)
    @test length(full_equations(simplify_busmodel(bm))) == 2
    bus = Bus(bm)
    toi = bus_on_slack(bus)
    toi["active power"]["electric power of swing"] = VIndex(2,:swing₊Pel)
    toi["active power"]["electric power of load"] = VIndex(2,:pqload₊P)
    toi["reactive power"]["electric power of load"] = VIndex(2,:pqload₊Q)
    isinteractive() && plottoi(toi)
    @reftest "Swing_and_load" toi
end
