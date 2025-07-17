"""
Test systems module for OpPoDyn.jl

This module provides pre-configured test systems for testing OpPoDyn functionality.
Similar to NetworkDynamics' ComponentLibrary.jl, this provides reusable network models.

Usage:
    (isinteractive() && @__MODULE__()==Main ? includet : include)("testsystems.jl")
    nw = TestSystems.load_ieee9bus()
    # Now you can test initialization, powerflow, etc.
"""
module TestSystems

using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using Graphs

export load_ieee9bus

"""
    load_ieee9bus()

Load the IEEE 9-bus test system.

Returns an uninitialized Network object with:
- 3 generator buses (SauerPai machines with AVR and governors)
- 3 load buses (ConstantYLoad)
- 3 transmission buses
- 9 branches (6 lines + 3 transformers)

This network is ready for powerflow solving and initialization testing.
"""
function load_ieee9bus()
    # Generator Bus Model
    function GeneratorBus(; machine_p=(;), avr_p=(;), gov_p=(;))
        @named machine = SauerPaiMachine(;
            vf_input=true,
            τ_m_input=true,
            S_b=100,
            V_b=1,
            ω_b=2π*60,
            R_s=0.000125,
            T″_d0=0.01,
            T″_q0=0.01,
            machine_p... # unpack machine parameters
        )
        @named avr = AVRTypeI(vr_min=-5, vr_max=5,
            Ka=20, Ta=0.2,
            Kf=0.063, Tf=0.35,
            Ke=1, Te=0.314,
            E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662,
            tmeas_lag=false,
            avr_p... # unpack AVR parameters
        )
        @named gov = TGOV1(R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5,
            gov_p... # unpack governor parameters
        )
        # generate the "injector" as combination of multiple components
        injector = CompositeInjector([machine, avr, gov]; name=:generator)

        # generate the MTKBus (i.e. the MTK model containg the busbar and the injector)
        mtkbus = MTKBus(injector)
    end

    # Load Bus Model
    function ConstantYLoadBus()
        @named load = ConstantYLoad()
        MTKBus(load; name=:loadbus)
    end

    # Generator parameters from RTDS datasheet
    gen1p = (;X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969, X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64)
    gen2p = (;X_ls=0.08958, X_d=0.8958, X′_d=0.1198, X″_d=0.11, X_q=0.8645, X′_q=0.1969, X″_q=0.11, T′_d0=6.00, T′_q0=0.535, H= 6.40)
    gen3p = (;X_ls=0.13125, X_d=1.3125, X′_d=0.1813, X″_d=0.18, X_q=1.2578, X′_q=0.2500, X″_q=0.18, T′_d0=5.89, T′_q0=0.600, H= 3.01)

    # Instantiate MTK models
    mtkbus1 = GeneratorBus(; machine_p=gen1p)
    mtkbus2 = GeneratorBus(; machine_p=gen2p)
    mtkbus3 = GeneratorBus(; machine_p=gen3p)
    mtkbus4 = MTKBus()
    mtkbus5 = ConstantYLoadBus()
    mtkbus6 = ConstantYLoadBus()
    mtkbus7 = MTKBus()
    mtkbus8 = ConstantYLoadBus()
    mtkbus9 = MTKBus()

    # Build NetworkDynamics components with powerflow models
    @named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.04))
    @named bus2 = Bus(mtkbus2; vidx=2, pf=pfPV(V=1.025, P=1.63))
    @named bus3 = Bus(mtkbus3; vidx=3, pf=pfPV(V=1.025, P=0.85))
    @named bus4 = Bus(mtkbus4; vidx=4)
    @named bus5 = Bus(mtkbus5; vidx=5, pf=pfPQ(P=-1.25, Q=-0.5))
    @named bus6 = Bus(mtkbus6; vidx=6, pf=pfPQ(P=-0.9, Q=-0.3))
    @named bus7 = Bus(mtkbus7; vidx=7)
    @named bus8 = Bus(mtkbus8; vidx=8, pf=pfPQ(P=-1.0, Q=-0.35))
    @named bus9 = Bus(mtkbus9; vidx=9)

    # Add initialization formulas for load buses
    vset_formula = @initformula :load₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
    add_initformula!(bus5, vset_formula)
    add_initformula!(bus6, vset_formula)
    add_initformula!(bus8, vset_formula)

    # Branch helper functions
    function piline(; R, X, B)
        @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
        MTKLine(pibranch)
    end
    function transformer(; R, X)
        @named xfmr = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
        MTKLine(xfmr)
    end

    # Define branches
    @named l45 = Line(piline(; R=0.0100, X=0.0850, B=0.1760), src=4, dst=5)
    @named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6)
    @named l57 = Line(piline(; R=0.0320, X=0.1610, B=0.3060), src=5, dst=7)
    @named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
    @named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
    @named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
    @named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
    @named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
    @named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)

    # Build the network
    vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9]
    edgefs = [l45, l46, l57, l69, l78, l89, t14, t27, t39]

    # Return uninitialized network
    return Network(vertexfs, edgefs)
end

end # module TestSystems