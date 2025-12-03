using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit

function Generator(; machine_p=(;), name)
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
    @named avr = AVRTypeI(; vr_min=-5, vr_max=5,
        Ka=20, Ta=0.2,
        Kf=0.063, Tf=0.35,
        Ke=1, Te=0.314,
        E1=3.3, Se1=0.6602, E2=4.5, Se2=4.2662,
        tmeas_lag=false,
    )
    @named gov = TGOV1(; R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)
    ## generate the "injector" as combination of multiple components
    injector = CompositeInjector([machine, avr, gov]; name)
end

gen1p = (;X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969, X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64)
gen2p = (;X_ls=0.08958, X_d=0.8958, X′_d=0.1198, X″_d=0.11, X_q=0.8645, X′_q=0.1969, X″_q=0.11, T′_d0=6.00, T′_q0=0.535, H= 6.40)
@named gen1 = Generator(machine_p=gen1p)
@named gen2 = Generator(machine_p=gen1p)
@named pqload = Library.PQConstraint(P=-1.0, Q=-0.5)
@named shunt = Library.ConstantYLoad(G=0.1, B=0.2)

@named busmod = MTKBus([gen1, gen2, pqload, shunt])
@named main = compile_bus(busmod, vidx=1, pf=pfSlack(V=1.0, δ=0.0))

line = compile_line(MTKLine(PiLine(R=0.01, name=:piline)), src=1, dst=2)
load = pfPQ(P=-1, Q=-0.1, vidx=2)

nw = Network([main, load], line)
s0 = initialize_from_pf!(nw; subverbose=true)
break

dump_initial_state(main)

multipf = Dict(
    :gen1 => SlackType(V=1.0),
    :gen2 => PVType(P=0.8, V=1.0),
    :pqload => PQType(P=-1.0, Q=-0.5),
    :shunt => YType(G=0.1, B=0.2),
)
set_pfmodel!(main, multipf)

multipf

to_vertexmodel
vm = main

