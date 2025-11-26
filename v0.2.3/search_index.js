var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "General",
    "title": "General",
    "category": "page",
    "text": "(Image: Build Status) (Image: Chat on Slack.) (Image: Get your Slack invitation.) (Image: Code on Github.)"
},

{
    "location": "index.html#PowerDynamics.jl-Dynamic-Power-System-Analysis-in-Julia-1",
    "page": "General",
    "title": "PowerDynamics.jl - Dynamic Power System Analysis in Julia",
    "category": "section",
    "text": "This package provides all the tools you need to create a dynamic power grid model and analyze it.The source code is licensed under GPLv3 and published on github.These Docs have been built with the following version of the sub-packages:using Pkg\nenv = Pkg.Types.EnvCache()\nsubpackages = copy(env.project[\"deps\"])\npop!(subpackages, \"Reexport\") # ignore Reexport here\nif haskey(subpackages, \"Documenter\")\n  pop!(subpackages, \"Documenter\") # ignore Documenter here\nend\nuuids = map(Base.UUID, values(subpackages))\nfunction printversions()\n  for uuid in uuids\n    pkg = Pkg.Types.manifest_info(env, uuid)\n    println(rpad(\" * $(pkg[\"name\"]) \", 30, \".\"), \" $(pkg[\"version\"])\")\n  end\nendprintversions() # hide"
},

{
    "location": "index.html#Installation-1",
    "page": "General",
    "title": "Installation",
    "category": "section",
    "text": "The installation can be done via the new package manager. Either use]add PowerDynamicsor copyusing Pkg; Pkg.add(\"PowerDynamics\")Please note that PowerDynamics.jl is a fast developing library whose API is not settled yet. In order to ensure that your old code will still work in the future while using the latest version of PowerDynamics.jl for your new code, we strongly recommend the usage of environments. Please check out this video from the introduction of Pkg3, where environments are introduced, too."
},

{
    "location": "index.html#Compatibility-1",
    "page": "General",
    "title": "Compatibility",
    "category": "section",
    "text": "PowerDynamics.jl is written for Julia 1.0 and above. We will quickly switch to new Julia version as they come out, but support older versions and enable long transition periods for users. Julia versions 0.x are not supported."
},

{
    "location": "index.html#Usage-1",
    "page": "General",
    "title": "Usage",
    "category": "section",
    "text": "Generally, we distinguish three types of user for PowerDynamics.jl:Grid Modeler\nGrid Component Developer\nPowerDynamics.jl Developer"
},

{
    "location": "index.html#Grid-Modeler-1",
    "page": "General",
    "title": "Grid Modeler",
    "category": "section",
    "text": "Your Goal is to use PowerDynamics.jl to model your grid of preference. You don\'t want to implement new types of nodes.We recommend you to choose your favorite example from PowerDynamicsExamples, read Node Types and try to understand it. That should give you the kickstart you need. If you have any questions, contact us."
},

{
    "location": "index.html#Grid-Component-Developer-1",
    "page": "General",
    "title": "Grid Component Developer",
    "category": "section",
    "text": "Your Goal is to use PowerDynamics.jl to develop types of nodes, e.g. new control schemes for inverters or new descriptions of synchronous machines.After going through the introduction for a Grid Modeler, we recommend that you read through Dynamics Types and Custom Node Types and try to implement a new node type for an example grid. With that, you should have all the tools you need. If you have any questions, contact us."
},

{
    "location": "index.html#PowerDynamics.jl-Developer-1",
    "page": "General",
    "title": "PowerDynamics.jl Developer",
    "category": "section",
    "text": "Your Goal is to extend PowerDynamics.jl with new fundamental functionalities.After going throught the introduction for a Grid Modeler and a Grid Component Developer, read through the code where hopefully all of this documentation will helpful for you. Afterwards, it\'s probably best to open an issue explainng the idea you want to implement and we can discuss how you can transform your idea into a pull request."
},

{
    "location": "language_conventions.html#",
    "page": "Language & Conventions",
    "title": "Language & Conventions",
    "category": "page",
    "text": ""
},

{
    "location": "language_conventions.html#Language-and-Conventions-1",
    "page": "Language & Conventions",
    "title": "Language & Conventions",
    "category": "section",
    "text": "Generally, variables are miniscule (e.g. u, i, ω) and parameters are capitalized (e.g. H, D, P, Ω). As it is common to use greek letters for modeling equations and Julia supports Unicode, greek letters are used within the Code, e.g. Ω and ω in PowerDynBase.SwingEq. If you don\'t want to use greek keyboard (which I am currently switching to) you can simply type the latex representating \\Omega and Atom can complete it with Ω using Tab."
},

{
    "location": "language_conventions.html#List-of-symbols-and-corresponding-names-1",
    "page": "Language & Conventions",
    "title": "List of symbols and corresponding names",
    "category": "section",
    "text": "Symbol (Code) Symbol (Math) Name within PowerDynamics.jl Common alternative names\n  node bus, vertex\n  grid network, power grid, power network\n y_ab = y_ba admittance between nodes a and b \nLY Y^L admittance laplacian (nodal) admittance matrix\nt t time \nim j imaginary element sqrt-1\nu = v \\cdot exp(im*φ) u = v cdot e^jφ complex voltage \nv v voltage magnitude absolute voltage\nφ phi voltage angle \ni_c = i \\cdot exp(im*δ) i_c = i cdot e^jdelta nodal complex current \ni i magnitude of the current \nδ delta angle of the current \ns = p + im*q s = p + jq complex power \np p real power active power\nq q imaginary power reactive power"
},

{
    "location": "language_conventions.html#List-of-modeling-conventions-1",
    "page": "Language & Conventions",
    "title": "List of modeling conventions",
    "category": "section",
    "text": "Counting of nodes starts at 1.\nRanges of nodes are mathematical, i.e. they include the first and the last element. For example sum_k=3^6 sums over 3, 4, 5, and 6.\nFor now, no selfadmittance is allowed, i.e. y_aa = 0 for all nodes a.\nThe admittance laplacian uses the following definition (convention from wikipedia)Y^L_ab = begincases\n  sum_c y_ac  textif  a=b \n  -y_ab  textotherwise\nendcasesThe nodal complex current is calculated asi_c_a = sum_b LY_ab u_b The complex power is calculated as (with ^* as complex comjucation)s_a = u_a cdot i_c_a^*"
},

{
    "location": "node_dynamics_types.html#",
    "page": "Dynamics Types",
    "title": "Dynamics Types",
    "category": "page",
    "text": ""
},

{
    "location": "node_dynamics_types.html#Dynamics-Types-1",
    "page": "Dynamics Types",
    "title": "Dynamics Types",
    "category": "section",
    "text": ""
},

{
    "location": "node_dynamics_types.html#PowerDynBase.AbstractNodeDynamics",
    "page": "Dynamics Types",
    "title": "PowerDynBase.AbstractNodeDynamics",
    "category": "type",
    "text": "Abstract super type for all abstract node dynamics types.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.OrdinaryNodeDynamics",
    "page": "Dynamics Types",
    "title": "PowerDynBase.OrdinaryNodeDynamics",
    "category": "type",
    "text": "OrdinaryNodeDynamics(;rhs, n_int)\n\nThe type representing the dynamics of a node that is described via ODEs.\n\nEach node a has the complex voltage u and n real internal variables y_1 dots y_n, so it generally describes a system of ordinary differential equation as\n\nfracdu_adt = f_u(u_a i_c_a y_1 dots y_n) \nfracdy_akdt = f_k(u_a i_c_a y_1 dots y_n)quad forall k = 1 dots n\n\nf is represented by rhs field of OrdinaryNodeDynamics.\n\nthe general signature of rhs is\n\nrhs(dint_dt::AbstractVector,\n    u::Complex,\n    i::Complex,\n    int::AbstractVector,\n    t,\n    )::Complex\n\nInput\nu is the complex voltage u\ni is the complex current i\nint is the array of internal variables y_1 dots y_n\nt is the time t\nOutput\nthe (complex) return value describes fracdudt\nrhs writes values in dint_dt describing the left-hand side fracdy_1dt dots fracdy_ndt\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.OrdinaryNodeDynamicsWithMass",
    "page": "Dynamics Types",
    "title": "PowerDynBase.OrdinaryNodeDynamicsWithMass",
    "category": "type",
    "text": "OrdinaryNodeDynamicsWithMass(;rhs, n_int, m_u, m_int)\n\nThe type representing the dynamics of a node that is described via ODEs.\n\nEach node a has the complex voltage u and n (= n_int) real internal variables y_1 dots y_n, so it generally describes a system of ordinary differential equation with a voltage mass m_u and internal masses m^int_1 dots m^int_n as\n\nm_ufracdu_adt = f_u(u_a i_c_a y_1 dots y_n) \nm^int_kfracdy_akdt = f_k(u_a i_c_a y_1 dots y_n)quad forall k = 1 dots n\n\nAs we assume that all masses are binary (either 1, or 0), that means, one can implement semi-explicit differential algebraic equations with this node dynamics type. f is represented by rhs field of OrdinaryNodeDynamics.\n\nthe general signature of rhs is\n\nrhs(dint_dt::AbstractVector,\n    u::Complex,\n    i::Complex,\n    int::AbstractVector,\n    t,\n    )::Complex\n\nInput\nu is the complex voltage u\ni is the complex current i\nint is the array of internal variables y_1 dots y_n\nt is the time t\nOutput\nthe (complex) return value describes fracdudt\nrhs writes values in dint_dt describing the left-hand side fracdy_1dt dots fracdy_ndt\n\nThe binary masses are:\n\nm_u is the boolean value for m_u\nm_int is the array of boolean values for m^int_1 dots m^int_n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Node-Dynamics-1",
    "page": "Dynamics Types",
    "title": "Node Dynamics",
    "category": "section",
    "text": "In this section, the implemented general types of node dynamics (which are all subtypes of PowerDynBase.AbstractNodeDynamics) with the corresponding helper functions and constants are introduced. The documentation is done for each type below. The main types are:PowerDynBase.OrdinaryNodeDynamics\nPowerDynBase.OrdinaryNodeDynamicsWithMassPowerDynBase.AbstractNodeDynamics\nPowerDynBase.OrdinaryNodeDynamics\nPowerDynBase.OrdinaryNodeDynamicsWithMass"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.internalsymbolsof",
    "page": "Dynamics Types",
    "title": "PowerDynBase.internalsymbolsof",
    "category": "function",
    "text": "Get the symbols representing the internal variables of the node.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Helper-Functions-1",
    "page": "Dynamics Types",
    "title": "Helper Functions",
    "category": "section",
    "text": "to be donePowerDynBase.internalsymbolsof"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.GridDynamics",
    "page": "Dynamics Types",
    "title": "PowerDynBase.GridDynamics",
    "category": "type",
    "text": "Abstract super type for all abstract grid dynamics types.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.OrdinaryGridDynamics",
    "page": "Dynamics Types",
    "title": "PowerDynBase.OrdinaryGridDynamics",
    "category": "type",
    "text": "struct OrdinaryGridDynamics <: AbstractOrdinaryGridDynamics\n    rhs::NetworkRHS\nend\n\nThe data structure that contains all the information necessary for a power grid model that can be described as an ordinary differential equation. In this case, only the PowerDynBase.NetworkRHS is necessary.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.OrdinaryGridDynamicsWithMass",
    "page": "Dynamics Types",
    "title": "PowerDynBase.OrdinaryGridDynamicsWithMass",
    "category": "type",
    "text": "struct OrdinaryGridDynamicsWithMass <: AbstractAlgebraicGridDynamics\n    rhs::NetworkRHS\n    masses::AbstractVector{Bool} # diagonal part of the mass matrix, off-diagonal is assumed to be 0 anyway\nend\n\nThe data structure that contains all the information necessary for a power grid model that can be described as an ordinary differential equation with masses, i.e. a semi-explicit differential algebraic equation. rhs is the PowerDynBase.NetworkRHS. masses is a 1-dimensional array representing the diagonal entries of the mass matrix. The off-diagonal entries are assumed to be 0. masses can only contain boolean values representing: true the equation is treated as a ordinary differential eqation and false the equation is treated as an algebraic constraint on the state variables.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.AlgebraicGridDynamics",
    "page": "Dynamics Types",
    "title": "PowerDynBase.AlgebraicGridDynamics",
    "category": "type",
    "text": "struct AlgebraicGridDynamics <: AbstractAlgebraicGridDynamics\n    rhs::NetworkRHS\n    differentials::AbstractVector{Bool} # boolean values whether there a variable is a differential\nend\n\nThe data structure that contains all the information necessary for a power grid model that can be described as an differential algebraic equation. rhs is the PowerDynBase.NetworkRHS. differentials is a 1-dimensional array of boolean values. A true entry means the corresponding variable is dynamic and has a derivative variable. A false entry means the corresponding variable is defined by an algebraic constraint only.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Grid-Dynamics-1",
    "page": "Dynamics Types",
    "title": "Grid Dynamics",
    "category": "section",
    "text": "Analogously, for each of the node dynamics types exists a corresponding grid dynamics type, that represents the dynamics of a whole power grid model. They are all subtypes of PowerDynBase.GridDynamics These are:PowerDynBase.OrdinaryGridDynamics\nPowerDynBase.OrdinaryGridDynamicsWithMass\nPowerDynBase.AlgebraicGridDynamicsPowerDynBase.GridDynamics\nPowerDynBase.OrdinaryGridDynamics\nPowerDynBase.OrdinaryGridDynamicsWithMass\nPowerDynBase.AlgebraicGridDynamics"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.NetworkRHS",
    "page": "Dynamics Types",
    "title": "PowerDynBase.NetworkRHS",
    "category": "type",
    "text": "struct NetworkRHS{T, M} <: AbstractNetworkFunction{T, M}\n    nodes::AbstractVector{T}\n    LY::M\n    numnodes\n    systemsize\n    intrange # unitrange telling me where I find the internal dynamic variables\n    nodalintranges # unit ranges to find the internal variables for each node in the full length of internal variables\nend\n\nRepresenting the full dynamics of the power grid.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#NetworkRHS-1",
    "page": "Dynamics Types",
    "title": "NetworkRHS",
    "category": "section",
    "text": "The logic of building a full power grid model (i.e. a subtype of PowerDynBase.GridDynamics) is encoded in PowerDynBase.NetworkRHS. Currently, the docs are a bit thin here.PowerDynBase.NetworkRHS"
},

{
    "location": "node_dynamics_types.html#Helper-Functions-2",
    "page": "Dynamics Types",
    "title": "Helper Functions",
    "category": "section",
    "text": "to be done"
},

{
    "location": "node_types.html#",
    "page": "Node Types",
    "title": "Node Types",
    "category": "page",
    "text": ""
},

{
    "location": "node_types.html#PowerDynBase.AbstractNodeParameters",
    "page": "Node Types",
    "title": "PowerDynBase.AbstractNodeParameters",
    "category": "type",
    "text": "Abstract super type for all node parameter types.\n\n\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.PQAlgebraic",
    "page": "Node Types",
    "title": "PowerDynBase.PQAlgebraic",
    "category": "type",
    "text": "PQAlgebraic(;S)\n\nA node type that locally fixes the active (P) and reactive power (Q) output of the node.\n\nKeyword Arguments\n\nS = P + Q*im: the complex power output\n\nMathematical Representation\n\nUsing PQAlgebraic for node a applies the equation\n\n0 = S_a - u_a cdot i_a^*\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.PVAlgebraic",
    "page": "Node Types",
    "title": "PowerDynBase.PVAlgebraic",
    "category": "type",
    "text": "PVAlgebraic(;P,V)\n\nA node type that locally fixes the active power (P) and the voltage magnitude (V) of the node.\n\nKeyword Arguments\n\nP: the active (real) power output\nV: voltage magnitude\n\nMathematical Representation\n\nUsing PVAlgebraic for node a applies the equations\n\n0 = P_a - Releft(u_a cdot i_a^*right) \n0 = V_a - leftu_aright\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.SlackAlgebraic",
    "page": "Node Types",
    "title": "PowerDynBase.SlackAlgebraic",
    "category": "type",
    "text": "SlackAlgebraic(;U)\n\nA node type that locally fixes the complex voltage (U) of the node.\n\nAs the complex voltage can be represented as U=Ve^iphi, this is equivlant to fixing the voltage magnitude V and the angle phi.\n\nKeyword Arguments\n\nU: the complex voltage\n\nMathematical Representation\n\nUsing SlackAlgebraic for node a applies the equation\n\n0 = U_a - u_a\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.SwingEq",
    "page": "Node Types",
    "title": "PowerDynBase.SwingEq",
    "category": "type",
    "text": "SwingEq(;H, P, D, Ω)\n\nA node type that applies the swing equation to the frequency/angle dynamics and keeps the voltage magnitude as is.\n\nAdditionally to u, it has the internal dynamic variable omega representing the frequency of the rotator relative to the grid frequency Omega, i.e. the real frequency omega_r of the rotator is given as omega_r = Omega + omega.\n\nKeyword Arguments\n\nH: inertia\nP: active (real) power output\nD: damping coefficient\nΩ: rated frequency of the power grid, often 50Hz\n\nMathematical Representation\n\nUsing SwingEq for node a applies the equations\n\nfracdu_adt = i u_a  omega_a \nfracH2piOmegafracdomega_adt = P_a - D_aomega_a - Releft(u_a cdot i_a^*right)\n\nwhich is equivalent to\n\nfracdphi_adt = omega \nv = v(t=0) = textconst \nfracH2piOmegafracdomega_adt = P_a - D_aomega_a - Releft(u_a cdot i_a^*right)\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.SwingEqLVS",
    "page": "Node Types",
    "title": "PowerDynBase.SwingEqLVS",
    "category": "type",
    "text": "SwingEqLVS(;H, P, D, Ω, Γ, V)\n\nA node type that applies the swing equation to the frequency/angle dynamics and has a linear voltage stability (LVS) term.\n\nAdditionally to u, it has the internal dynamic variable omega representing the frequency of the rotator relative to the grid frequency Omega, i.e. the real frequency omega_r of the rotator is given as omega_r = Omega + omega.\n\nKeyword Arguments\n\nH: inertia\nP: active (real) power output\nD: damping coefficient\nΩ: rated frequency of the power grid, often 50Hz\nΓ: voltage stability coefficient\nV: set voltage, usually 1\n\nMathematical Representation\n\nUsing SwingEq for node a applies the equations\n\nfracdu_adt = i u_a omega - fracuu Γ_a  (v_a - V_a) \nfracH2piOmegafracdomega_adt = P_a - D_aomega_a - Releft(u_a cdot i_a^*right)\n\nwhich is equivalent to\n\nfracdphi_adt = omega_a \nfracdv_adt = - Γ_a  (v_a - V_a) \nfracH2piOmegafracdomega_adt = P_a - D_aomega_a - Releft(u_a cdot i_a^*right)\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.FourthEq",
    "page": "Node Types",
    "title": "PowerDynBase.FourthEq",
    "category": "type",
    "text": "FourthEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q)\n\nA node type that applies the 4th-order synchronous machine model with frequency/angle and voltage dynamics.\n\nAdditionally to u, it has the internal dynamic variables\n\nomega representing the frequency of the rotator relative to the grid frequency Omega, i.e. the real frequency omega_r of the rotator is given as omega_r = Omega + omega and\ntheta representing the relative angle of the rotor with respect to the voltage angle phi.\n\nKeyword Arguments\n\nH: inertia\nP: active (real) power output\nD: damping coefficient\nΩ: rated frequency of the power grid, often 50Hz\nT_d_dash: time constant of d-axis\nT_q_dash: time constant of q-axis\nX_d_dash: transient reactance of d-axis\nX_q_dash: transient reactance of q-axis\nX_d: reactance of d-axis\nX_d: reactance of q-axis\n\nMathematical Representation\n\nUsing FourthEq for node a applies the equations\n\n    u = -je_c e^jtheta = -j(e_d + je_q)e^jtheta\n    e_c= e_d + je_q = jue^-jtheta\n    i  = -jie^jtheta = -j(i_d+ j i_q )e^jtheta = Y^L cdot u \n    i_c= i_d + ji_q = jie^-jtheta\n    p = Re (i^* u)\n\nThe fourth-order equations read (according to Sauer, p. 140, eqs. (6110)-(6114)) and p. 35 eqs(3.90)-(3.91)\n\n    fracdthetadt = omega \n     fracdomegadt = P-Domega - p -(x_q-x_d)i_d i_q\n    fracd e_qdt = frac1T_d (- e_q - (x_d - x_d) i_d+ e_f) \n    fracd e_ddt = frac1T_q (- e_d + (x_q - x_q) i_q)  \n\nWith the PowerDynamics.jl \\time{naming conventions} of i and u they read as\n\n   dot u = fracddt(-j e_c e^jtheta)=-j(dot e_d + jdot e_q)e^jtheta + ujomega\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.VSIMinimal",
    "page": "Node Types",
    "title": "PowerDynBase.VSIMinimal",
    "category": "type",
    "text": "VSIMinimal(;τ_P,τ_Q,K_P,K_Q,E_r,P,Q)\n\nA node type that applies the frequency and voltage droop control to control the frequency and voltage dynamics.\n\nAdditionally to u, it has the internal dynamic variable omega representing the frequency of the rotator relative to the grid frequency Omega, i.e. the real frequency omega_r of the rotator is given as omega_r = Omega + omega.\n\nKeyword Arguments\n\nτ_p: time constant active power measurement\nτ_Q: time constant reactive power measurement\nK_P: droop constant frequency droop\nK_Q: droop constant voltage droop\nV_r: reference/ desired voltage\nP: active (real) power infeed\nQ: reactive (imag) power infeed\n\nMathematical Representation\n\nUsing VSIMinimal for node a applies the equations\n\ndotphi_a=omega_a\n dotomega_a=frac1tau_Pa-omega_a-K_Pa (Releft(u_a cdot i_a^*right)-P_refa)\ntau_Qdotv_a=-v_a+V_ref-K_Qa (Imleft(u_a cdot i_a^*right)-Q_refa)\n dotu_a=dotv_ae^jphi+jomega_a u_a\n\n```\n\n\n\n"
},

{
    "location": "node_types.html#PowerDynBase.VSIVoltagePT1",
    "page": "Node Types",
    "title": "PowerDynBase.VSIVoltagePT1",
    "category": "type",
    "text": "VSIVoltagePT1(;τ_v,τ_P,τ_Q,K_P,K_Q,E_r,P,Q)\n\nA node type that applies the frequency and voltage droop control to control the frequency and voltage dynamics.\n\nAdditionally to u, it has the internal dynamic variable omega representing the frequency of the rotator relative to the grid frequency Omega, i.e. the real frequency omega_r of the rotator is given as omega_r = Omega + omega.\n\nKeyword Arguments\n\nτ_v: time constant voltage control delay\nτ_p: time constant active power measurement\nτ_Q: time constant reactive power measurement\nK_P: droop constant frequency droop\nK_Q: droop constant voltage droop\nV_r: reference/ desired voltage\nP: active (real) power infeed\nQ: reactive (imag) power infeed\n\nMathematical Representation\n\nUsing VSIVoltagePT1 for node a applies the equations\n\ndotphi_a=omega_a\n dotomega_a=frac1tau_Pa-omega_a-K_Pa (Releft(u_a cdot i_a^*right)-P_refa)\n tau_vdotv_a=-v_a+V_ref-K_Qa(q_ma-Q_refa)\n tau_Q dotq_ma=-q_ma+Imleft(u_a cdot i_a^*right)\n dotu_a=dotv_ae^jphi+jomega_a u_a\n\n\n\n"
},

{
    "location": "node_types.html#Node-Types-1",
    "page": "Node Types",
    "title": "Node Types",
    "category": "section",
    "text": "The currently implementes node types arePurely Algebraic:\nPowerDynBase.PQAlgebraic (PQ-bus)\nPowerDynBase.PVAlgebraic (PV-bus)\nPowerDynBase.SlackAlgebraic (Slack-bus / Vφ-bus)\nSynchronous Machine Models:\nPowerDynBase.SwingEq (2nd order)\nPowerDynBase.SwingEqLVS (2nd order with an additional term for numerical voltage stability)\nPowerDynBase.FourthEq (4th order)\nVoltage Source Inverters:\nPowerDynBase.VSIMinimal\nPowerDynBase.VSIVoltagePT1They are all subtypes of PowerDynBase.AbstractNodeParameters.PowerDynBase.AbstractNodeParameters\nPQAlgebraic\nPVAlgebraic\nSlackAlgebraic\nSwingEq\nSwingEqLVS\nFourthEq\nVSIMinimal\nVSIVoltagePT1"
},

{
    "location": "custom_node_types.html#",
    "page": "Custom Node Types",
    "title": "Custom Node Types",
    "category": "page",
    "text": ""
},

{
    "location": "custom_node_types.html#PowerDynBase.@DynamicNode",
    "page": "Custom Node Types",
    "title": "PowerDynBase.@DynamicNode",
    "category": "macro",
    "text": "Macro for creating a new type of dynmic nodes.\n\nSyntax Description:\n\n@DynamicNode MyNewNodeName(Par1, Par2, ...) <: NodeDynamicsType(N1, N2, ...) begin\n    [all prepratory things that need to be run just once]\nend [[x1, dx1], [x2, dx2]] begin\n    [the actual dynamics equation]\n    [important to set the output variables]\nend\n\nwhere MyNewNodeName is the name of the new type of dynamic node, Par1, Par2, ... are the names of the parameters, NodeDynamicsType the the node dynamics type (e.g. OrdinaryNodeDynamics or OrdinaryNodeDynamicsWithMass), N1, N1, ... the parameters of the dynamics type, x1, x2, ... the internal variables of the node and dx1, dx2, ... the corresponding differentials.\n\nIn the first block, the preparation code that needs to be run only once is inserted. Finally, the second block contains the dynamics description, where it\'s important that the output variables need to be set. In case of OrdinaryNodeDynamics and OrdinaryNodeDynamicsWithMass, these are du and the differentials of the internal variables (here dx1, dx2).\n\nBelow are two examples:\n\n@DynamicNode SwingEqParameters(H, P, D, Ω) <: OrdinaryNodeDynamics() begin\n    @assert D > 0 \"damping (D) should be >0\"\n    @assert H > 0 \"inertia (H) should be >0\"\n    Ω_H = Ω * 2pi / H\nend [[ω, dω]] begin\n    p = real(u * conj(i_c))\n    dϕ = ω # dϕ is only a temp variable that Julia should optimize out\n    du = u * im * dϕ\n    dω = (P - D*ω - p)*Ω_H\nend\n\n@DynamicNode SlackAlgebraicParameters(U) <: OrdinaryNodeDynamicsWithMass(m_u=false, m_int=no_internal_masses) begin\nend [] begin\n        du = u - U\nend\n\n\n\n\n\n"
},

{
    "location": "custom_node_types.html#Custom-Node-Types-1",
    "page": "Custom Node Types",
    "title": "Custom Node Types",
    "category": "section",
    "text": "To define your own Node Types, use the PowerDynBase.@DynamicNode macro. The new node type will be a subtype of PowerDynBase.AbstractNodeParameters.@DynamicNode"
},

{
    "location": "states_solutions.html#",
    "page": "States, Solutions & Plotting",
    "title": "States, Solutions & Plotting",
    "category": "page",
    "text": ""
},

{
    "location": "states_solutions.html#States,-Solutions-and-Plotting-1",
    "page": "States, Solutions & Plotting",
    "title": "States, Solutions & Plotting",
    "category": "section",
    "text": "In order to properly interact with the state space of the power grid model, we defined two data structures [PowerDynBase.State] and [PowerDynSolve.GridSolution]."
},

{
    "location": "states_solutions.html#PowerDynBase.State",
    "page": "States, Solutions & Plotting",
    "title": "PowerDynBase.State",
    "category": "type",
    "text": "\n    State(base; t=nothing)\n    State(grid, vec; t=nothing)\n\n\nEncode the information on the value of a state vector at a particular time point.\n\nKeyword Arguments\n\nbase is an instance of a BaseState, essentially it contains the state   vector and the complete rhs of the system. Instead of base, you can also   directly use a GridDynamics instance grid and a properly sized   state vector vec to instantiate a State.\nt is a time point associated to the base. It defaults to nothing.\n\nIndexing\n\nConcerning the indexing, a State object s basically behaves like a an array. There are plenty of convenient ways to access its contents at a node j by using a particular symbol:\n\ns[j, :u]: complex voltage\ns[j, :v]: voltage magnitude\ns[j, :φ]: voltage angle\ns[j, :i]: complex nodal current\ns[j, :iabs]: nodal current magnitude\ns[j, :δ]: nodal current angle\ns[j, :s]: apparent power\ns[j, :p]: real power\ns[j, :q]: reactive power\n\nCurrently, setting the state value is only implemented for u and v, the other quantities are derived automatically.\n\nWhen a node j has internal variables, you can access (and set) the k-th internal variable by calling\n\ns[j, :int, k].\n\nThe internal variables can be also directly accessed with symbols, i.e.\n\ns[j, :ω]\n\nreturns the frequency ω at node j. To find out the proper symbol, the easiest way is to look into the docs of the corresponding AbstractNodeParameters subtype, check the output of internalsymbolsof or simply look at the output of println:\n\njulia> internalsymbolsof(SwingEq(H=2, P=3, D=4, Ω=5))\n1-element Array{Symbol,1}:\n :ω\n\njulia> println(SwingEq(H=2, P=3, D=4, Ω=5))\nSwingEq[:ω](H=2, P=3, D=4, Ω=5)\n\n\n\n"
},

{
    "location": "states_solutions.html#States-1",
    "page": "States, Solutions & Plotting",
    "title": "States",
    "category": "section",
    "text": "State"
},

{
    "location": "states_solutions.html#PowerDynSolve.GridSolution",
    "page": "States, Solutions & Plotting",
    "title": "PowerDynSolve.GridSolution",
    "category": "type",
    "text": "struct GridSolution <: AbstractGridSolution\n    dqsol::AbstractTimeseriesSolution\n    griddynamics::GridDynamics\nend\n\nThe data structure interfacing to the solution of the differntial equations of a power grid. Normally, it is not created by hand but return from PowerDynSolve.solve.\n\nAccessing the solution in a similar interface as PowerDynBase.State.\n\nFor some grid solution sol, one can access the variables as\n\nsol(t, n, s)\n\nwhere t is the time (either float or array), n the node number(s) (either integer, array, range (e.g. 2:3) or colon (:, for all nodes)), and s is the symbol represnting the chosen value. s can be either: :v, :φ, :i, :iabs, :δ, :s, :p, :q, or the symbol of the internal variable of the nodes. The meaning of the symbols derives from the conventions of PowerDynamics.jl.\n\nFinally, one can access the a-th internal variable of a node by using sol(t, n, :int, a).\n\nInterfacing the Plots.jl library via plotting recipes, that follow similar instructions as the direct access to the solution.\n\nFor some grid solution sol, one plot variables of the solution asin\n\nusing Plots\nplot(sol, n, s, plots_kwargs...)\n\nwhere n and s are as in the accessing of plots, and plots_kwargs are the keyword arguments for Plots.jl.\n\n\n\n\n\n"
},

{
    "location": "states_solutions.html#Solutions-and-Plotting-1",
    "page": "States, Solutions & Plotting",
    "title": "Solutions & Plotting",
    "category": "section",
    "text": "PowerDynSolve.GridSolution"
},

{
    "location": "error_types.html#",
    "page": "Error Types",
    "title": "Error Types",
    "category": "page",
    "text": ""
},

{
    "location": "error_types.html#PowerDynBase.PowerDynamicsError",
    "page": "Error Types",
    "title": "PowerDynBase.PowerDynamicsError",
    "category": "type",
    "text": "Abstract super type of all PowerDynamics.jl Errors.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynBase.NodeDynamicsError",
    "page": "Error Types",
    "title": "PowerDynBase.NodeDynamicsError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong during the node dynamics construction.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynBase.GridDynamicsError",
    "page": "Error Types",
    "title": "PowerDynBase.GridDynamicsError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong during the grid dynamics construction.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynBase.StateError",
    "page": "Error Types",
    "title": "PowerDynBase.StateError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong when creating or modifying states.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynSolve.GridSolutionError",
    "page": "Error Types",
    "title": "PowerDynSolve.GridSolutionError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong during when solving a power grid model.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#Error-Types-1",
    "page": "Error Types",
    "title": "Error Types",
    "category": "section",
    "text": "PowerDynBase.PowerDynamicsError\nPowerDynBase.NodeDynamicsError\nPowerDynBase.GridDynamicsError\nPowerDynBase.StateError\nPowerDynSolve.GridSolutionError"
},

{
    "location": "fullindex.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "fullindex.html#Index-1",
    "page": "Index",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "contact.html#",
    "page": "Contact",
    "title": "Contact",
    "category": "page",
    "text": "(Image: Build Status) (Image: Chat on Slack.) (Image: Get your Slack invitation.) (Image: Code on Github.)"
},

{
    "location": "contact.html#Contact-1",
    "page": "Contact",
    "title": "Contact",
    "category": "section",
    "text": "In case of questions, please submit an issue on github or ask on our slack channel (get your invitation here).If you don\'t want to contact us publicly, send an email to Tim Kittel (tim.kittel@elena-international.com) or Sabine Auer (sabine.auer@elena-international.com)."
},

]}
