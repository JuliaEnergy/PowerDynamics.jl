using OpPoDyn
using OpPoDyn.Library
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary: Blocks


systems = @named begin
    machine = DynawoMachine()
    excitation = Blocks.Constant(k=2.4659)
    pmech = Blocks.Constant(k=0.903)
    ωRef = Blocks.Constant(k=1.0)
end
eqs = [connect(machine.efdPu, excitation.output),
       connect(machine.PmPu, pmech.output),
       connect(machine.ωRefPu, ωRef.output)]
@named generator = ODESystem(eqs, t; systems)
bus = Bus((generator, generator.machine.terminal))
bus = pin_parameters(bus)
simp = structural_simplify_bus(bus)
unknowns(simp)



# complete(bus)
bus = pinparameters(bus, ModelingToolkit.defaults(bus))
# io = ([getproperty(bus, :busbar; namespace=false).i_r, getproperty(bus, :busbar; namespace=false).i_i], [])
io = (ModelingToolkit.unbound_inputs(bus), [])
simp = structural_simplify(bus, io)[1]

vert = ODEVertex(bus, [bus.busbar.i_r, bus.busbar.i_i], [bus.busbar.u_r, bus.busbar.u_i])
vert = ODEVertex(bus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i])

using OpPoDyn.Library

systems = @named begin
    machine = DynawoMachine()
    excitation = Blocks.Constant(k=2.4659)
    pmech = Blocks.Constant(k=0.903)
    ωRef = Blocks.Constant(k=1.0)
    trafo = DynawoFixedRatioTransformer()
    busbar = BusBar()
end
eqs = [connect(machine.efdPu, excitation.output),
       connect(machine.PmPu, pmech.output),
       connect(machine.ωRefPu, ωRef.output),
       connect(machine.terminal, trafo.terminal2),
       connect(trafo.terminal1, busbar.terminal)]
@named bus = ODESystem(eqs, t; systems)
bus = pin_parameters(bus)
simp = structural_simplify_bus(bus)
simp

structural_simplify(simp)

full_equations(simp)

ODEVertex(bus, )


# line model
@named branchA = DynawoPiLine(XPu=0.022522)
@named branchB = DynawoPiLine(XPu=0.04189)
@named line = Line(branchA, branchB)
line = pin_parameters(line)

# genbus model
@mtkmodel GenBus begin
    @components begin
        machine = DynawoMachine()
        excitation = Blocks.Constant(k=2.4659)
        pmech = Blocks.Constant(k=0.903)
        ωRef = Blocks.Constant(k=1.0)
        trafo = DynawoFixedRatioTransformer()
        busbar = BusBar()
    end
    @equations begin
        connect(machine.efdPu, excitation.output)
        connect(machine.PmPu, pmech.output)
        connect(machine.ωRefPu, ωRef.output)
        connect(machine.terminal, trafo.terminal2)
        connect(trafo.terminal1, busbar.terminal)
    end
end
@named bus = GenBus()
bus = pin_parameters(bus)

# slackbus model
# @named slack = SlackAlgebraic(;u_set_r=1)
@named slack = SlackDifferential(u_init_r=0.90081)


genf = vertex_model(bus)
slackf = vertex_model(slack)
linef = edge_model(line)

using Graphs
using NetworkDynamics
using OrdinaryDiffEq
g = path_graph(2)
nw = Network(g, [slackf, genf], linef)
u0 = NWState(nw)
u0.v[2, :machine₊θ]        = 0
u0.v[2, :machine₊λ_fPu]    = 0
u0.v[2, :machine₊λ_DPu]    = 0
u0.v[2, :machine₊λ_Q1Pu]   = 0
u0.v[2, :machine₊λ_Q2Pu]   = 0
u0.v[2, :machine₊ωPu]      = 1
u0.v[2, :machine₊idPu]     = 0
u0.v[2, :machine₊iqPu]     = 0
u0.v[2, :machine₊MdSat′Pu] = 0
u0.v[2, :machine₊ifPu]     = 0
u0.v[2, :machine₊iQ2Pu]    = 0
u0.v[2, :machine₊MqSat′Pu] = 0
u0.v[2, :machine₊iQ1Pu]    = 0

prob = ODEProblem(nw, uflat(u0), (0, 10), pflat(u0))
sol = solve(prob, Rodas5P())
sol = solve(prob, Rosenbrock23())


g = path_graph(2)
nw = Network(g, slackf, linef)
u0 = NWState(nw)

prob = ODEProblem(nw, uflat(u0), (0, 10), pflat(u0))
sol = solve(prob, Rodas5P())
