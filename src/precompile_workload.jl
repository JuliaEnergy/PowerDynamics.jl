# A 4-bus grid made of Library components only: a machine with AVR and governor on the slack
# bus, a junction, a load bus and a shunt bus. It is initialized from the powerflow on the sparse
# solver path, like a grid with more than 50 buses would be.
using .Library: SauerPaiMachine, AVRTypeI, TGOV1, VoltageDependentLoad, StaticShunt, PiLine

function _wl_run()
    @named machine = SauerPaiMachine(; vf_input=true, τ_m_input=true, R_s=0.000125, T″_d0=0.01,
        T″_q0=0.01, X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969,
        X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64)
    @named avr = AVRTypeI(; vr_min=-5, vr_max=5, Ka=20, Ta=0.2, Kf=0.063, Tf=0.35, Ke=1, Te=0.314,
        E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662, tmeas_lag=false)
    @named gov = TGOV1(; R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)
    gen = MTKBus(CompositeInjector([machine, avr, gov]; name=:generator))
    b1 = compile_bus(gen; vidx=1, pf=pfSlack(V=1.04), Vbase=16.5, name=:bus1)
    b2 = compile_bus(MTKBus(); vidx=2, Vbase=230, name=:bus2)
    load = VoltageDependentLoad(; name=:load, αP=1, αQ=2, Pset=-1.25, Qset=-0.5)
    b3 = compile_bus(MTKBus(load); vidx=3, pf=pfPQ(P=-1.25, Q=-0.5), Vbase=230, name=:bus3)
    shunt = StaticShunt(; name=:shunt, G=0, B=0.1)
    b4 = compile_bus(MTKBus(shunt); vidx=4, Vbase=230, name=:bus4)
    piline(; R, X, B) = MTKLine(PiLine(; R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, name=:piline))
    l12 = compile_line(piline(; R=0, X=0.0576, B=0); src=1, dst=2, name=:l12)
    l23 = compile_line(piline(; R=0.01, X=0.085, B=0.176); src=2, dst=3, name=:l23)
    l24 = compile_line(piline(; R=0.01, X=0.085, B=0.176); src=2, dst=4, name=:l24)
    nw = Network([b1, b2, b3, b4], [l12, l23, l24])
    initialize_from_pf(nw; verbose=false, sparsepf=true)
end
_wl_run()
