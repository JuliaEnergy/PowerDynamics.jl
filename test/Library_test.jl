using OpPoDyn
using OpPoDyn.Library
using OpPoDynTesting
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: ModelingToolkit as MTK
using OrderedCollections
using DiffEqCallbacks
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

@testset "Dynawo Machine test" begin
    # line model
    @named branchA = DynawoPiLine(XPu=0.022522)
    @named branchB = DynawoPiLine(XPu=0.04189)
    linem = LineModel(branchA, branchB)
    linef = Line(linem);

    # genbus model
    @mtkmodel GenBus begin
        @components begin
            machine = DynawoMachine()
            # machine = Swing(Pm_input=true)
            excitation = Blocks.Constant(k=2.4659)
            pmech = Blocks.Constant(k=0.903)
            ωRef = Blocks.Constant(k=1.0)
            trafo = DynawoFixedRatioTransformer()
            busbar = BusBar()
        end
        @equations begin
            connect(machine.efdPu, excitation.output)
            connect(machine.PmPu, pmech.output)
            # connect(machine.Pm, pmech.output)
            connect(machine.ωRefPu, ωRef.output)
            connect(machine.terminal, trafo.dst)
            connect(trafo.src, busbar.terminal)
        end
    end
    @named genbus = GenBus()
    # genbus = pin_parameters(genbus)
    genf = Bus(genbus; verbose=false);

    @named slack = SlackDifferential(u_init_r=0.90081)
    slackf = Bus(slack)

    g = path_graph(2)
    nw = Network(g, [slackf.compf, genf.compf], linef.compf)
    u0 = NWState(nw)
    u0.v[2, :machine₊θ]        = 1.2107
    u0.v[2, :machine₊λ_fPu]    = 1.1458
    u0.v[2, :machine₊λ_DPu]    = 0.89243
    u0.v[2, :machine₊λ_Q1Pu]   = -0.60044
    u0.v[2, :machine₊λ_Q2Pu]   = -0.60044
    u0.v[2, :machine₊ωPu]      = 1
    u0.v[2, :machine₊idPu]     = -0.91975
    u0.v[2, :machine₊iqPu]     = -0.39262
    u0.v[2, :machine₊MqSat′Pu] = 1.5292
    u0.v[2, :machine₊MdSat′Pu] = 1.5792
    u0.v[2, :machine₊ifPu]     = 1.4855
    u0.v[2, :machine₊iDPu]     = 0
    u0.v[2, :machine₊iQ2Pu]    = 0
    u0.v[2, :machine₊iQ1Pu]    = 0

    affect = function(int)
        p = NWParameter(int)
        p.v[2, :pmech₊k] = 0.923
        auto_dt_reset!(int)
        save_parameters!(int)
    end
    cb = PresetTimeCallback(10, affect)
    prob = ODEProblem(nw, uflat(u0), (0, 30), copy(pflat(u0)); callback=cb)
    sol = solve(prob, Rodas5P())

    plotspec = OrderedDict(
        "active power" => OrderedDict(
            "injection from bus" => VIndex(2, :busbar₊P),
            "PGenPu" => VIndex(2, :machine₊PGenPu)),
        "reactive power" => OrderedDict(
            "injection from bus" => VIndex(2, :busbar₊Q),
            "QGenPu" => VIndex(2, :machine₊QGenPu)),
        "voltage angle" => OrderedDict(
            "angle at bus" => VIndex(2, :busbar₊u_arg),
            "angle of rotor" => VIndex(2, :machine₊θ),
            "angle at slack" => VIndex(1, :busbar₊u_arg)),
        "voltage magnitude" => OrderedDict(
            "magnitude at bus" => VIndex(2, :busbar₊u_mag),
            "stator voltage" => VIndex(2, :machine₊UStatorPu),
            "magnitude at slack" => VIndex(1, :busbar₊u_mag)),
        "frequency" => OrderedDict(
            "frequency at machine" => VIndex(2, :machine₊ωPu)))
    toi = TrajectoriesOfInterest(sol, plotspec)
    isinteractive() && plottoi(toi)
    @reftest "dynawomachine_on_slack_prefstep" toi
end
