var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "PowerDynamics.jl",
    "title": "PowerDynamics.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#PowerDynamics.jl-Dynamic-Power-System-Analysis-in-Julia-1",
    "page": "PowerDynamics.jl",
    "title": "PowerDynamics.jl - Dynamic Power System Analysis in Julia",
    "category": "section",
    "text": "This package provides all the tools you need to create a dynamic power grid model and analyze it."
},

{
    "location": "index.html#Installation-1",
    "page": "PowerDynamics.jl",
    "title": "Installation",
    "category": "section",
    "text": "For now install all the packages using git directly. Please be aware of the order and make sure you have access rights. They will be registered with the official package managers upon publishing.TBD: add new installation instructions"
},

{
    "location": "index.html#Usage-1",
    "page": "PowerDynamics.jl",
    "title": "Usage",
    "category": "section",
    "text": "Generally, we distinguish three types of user for PowerDynamics.jl:Grid Modeler\nGrid Component Developer\nPowerDynamics.jl Developer"
},

{
    "location": "index.html#Grid-Modeler-1",
    "page": "PowerDynamics.jl",
    "title": "Grid Modeler",
    "category": "section",
    "text": "Your Goal is to use PowerDynamics.jl to model your grid of preference. You don\'t want to implement new types of nodes.We recommend you to choose your favorite example from PowerDynamicsExamples, read Node Types and try to understand it. That should give you the kickstart you need. If you have any questions, contact us."
},

{
    "location": "index.html#Grid-Component-Developer-1",
    "page": "PowerDynamics.jl",
    "title": "Grid Component Developer",
    "category": "section",
    "text": "Your Goal is to use PowerDynamics.jl to develop types of nodes, e.g. new control schemes for inverters or new descriptions of synchronous machines.After going through the introduction for a Grid Modeler, we recommend that you read through Node Dynamics Types and Custom Node Types and try to implement a new node type for an example grid. With that, you should have all the tools you need. If you have any questions, contact us."
},

{
    "location": "index.html#PowerDynamics.jl-Developer-1",
    "page": "PowerDynamics.jl",
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
    "page": "Node Dynamics Types",
    "title": "Node Dynamics Types",
    "category": "page",
    "text": ""
},

{
    "location": "node_dynamics_types.html#Node-Dynamics-Types-1",
    "page": "Node Dynamics Types",
    "title": "Node Dynamics Types",
    "category": "section",
    "text": "In this section, the implemented general types of node dynamics with the corresponding helper functions and constants are introduced. The documentation is done for each type below. The main types are:PowerDynBase.OrdinaryNodeDynamics\nPowerDynBase.OrdinaryNodeDynamicsWithMass"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.no_internal_differentials",
    "page": "Node Dynamics Types",
    "title": "PowerDynBase.no_internal_differentials",
    "category": "constant",
    "text": "A variable to be used when no internal differentials are present for a node dynamics type.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.no_internal_masses",
    "page": "Node Dynamics Types",
    "title": "PowerDynBase.no_internal_masses",
    "category": "constant",
    "text": "A variable to be used when no internal masses are present for a node dynamics type.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#PowerDynBase.OrdinaryNodeDynamicsWithMass",
    "page": "Node Dynamics Types",
    "title": "PowerDynBase.OrdinaryNodeDynamicsWithMass",
    "category": "type",
    "text": "OrdinaryNodeDynamicsWithMass(;rhs, n_int, m_u, m_int)\n\nThe type representing the dynamics of a node that is described via ODEs.\n\nEach node a has the complex voltage u and n (= n_int) real internal variables y_1 dots y_n, so it generally describes a system of ordinary differential equation with a voltage mass m_u and internal masses m^int_1 dots m^int_n as\n\nm_ufracdu_adt = f_u(u_a i_c_a y_1 dots y_n) \nm^int_kfracdy_akdt = f_k(u_a i_c_a y_1 dots y_n)quad forall k = 1 dots n\n\nAs we assume that all masses are binary (either 1, or 0), that means, one can implement semi-explicit differential algebraic equations with this node dynamics type.\n\nNo documentation found.\n\nMarkdown.@doc_str is a macro.\n\n# 1 method for macro \"@doc_str\":\n[1] @doc_str(__source__::LineNumberNode, __module__::Module, s::AbstractString, t...) in Markdown at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.0/Markdown/src/Markdown.jl:53\n\n\n\nThe binary masses are:\n\nm_u is the boolean value for m_u\nm_int is the array of boolean values for m^int_1 dots m^int_n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.convert-Tuple{Type{OrdinaryNodeDynamicsWithMass},OrdinaryNodeDynamics}",
    "page": "Node Dynamics Types",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{OrdinaryNodeDynamicsWithMass}, ::OrdinaryNodeDynamics)\n\nConversion of OrdinaryNodeDynamics to OrdinaryNodeDynamicsWithMass by assuming all masses are 1 (true).\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.convert-Tuple{Type{PowerDynBase.AlgebraicNodeDynamics},OrdinaryNodeDynamicsWithMass}",
    "page": "Node Dynamics Types",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{AlgebraicNodeDynamics}, ::OrdinaryNodeDynamicsWithMass)\n\nConversion of OrdinaryNodeDynamicsWithMass to AlgebraicNodeDynamics by transforming a RHS function into a root function.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.convert-Tuple{Type{PowerDynBase.AlgebraicNodeDynamics},OrdinaryNodeDynamics}",
    "page": "Node Dynamics Types",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{AlgebraicNodeDynamics}, ::OrdinaryNodeDynamics)\n\nConversion of OrdinaryNodeDynamics to AlgebraicNodeDynamics by going via OrdinaryNodeDynamicsWithMass.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.promote_rule-Tuple{Type{#s25} where #s25<:OrdinaryNodeDynamics,Type{#s24} where #s24<:OrdinaryNodeDynamicsWithMass}",
    "page": "Node Dynamics Types",
    "title": "Base.promote_rule",
    "category": "method",
    "text": "promote_rule(::Type{OrdinaryNodeDynamics}, ::Type{OrdinaryNodeDynamicsWithMass}) = OrdinaryNodeDynamicsWithMass\n\nOrdinaryNodeDynamics can be promoted to OrdinaryNodeDynamicsWithMass, see DPSABase.convert.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.promote_rule-Tuple{Type{OrdinaryNodeDynamicsWithMass},Type{PowerDynBase.AlgebraicNodeDynamics}}",
    "page": "Node Dynamics Types",
    "title": "Base.promote_rule",
    "category": "method",
    "text": "promote_rule(::Type{OrdinaryNodeDynamicsWithMass}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics\n\nOrdinaryNodeDynamicsWithMass can be promoted to AlgebraicNodeDynamics, see DPSABase.convert.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Base.promote_rule-Tuple{Type{OrdinaryNodeDynamics},Type{PowerDynBase.AlgebraicNodeDynamics}}",
    "page": "Node Dynamics Types",
    "title": "Base.promote_rule",
    "category": "method",
    "text": "promote_rule(::Type{OrdinaryNodeDynamics}, ::Type{AlgebraicNodeDynamics}) = AlgebraicNodeDynamics\n\nOrdinaryNodeDynamics can be promoted to AlgebraicNodeDynamics, see DPSABase.convert.\n\n\n\n\n\n"
},

{
    "location": "node_dynamics_types.html#Documentation-1",
    "page": "Node Dynamics Types",
    "title": "Documentation",
    "category": "section",
    "text": "Modules = [PowerDynBase]\nPages   = [\"NodeDynamicsBase.jl\"]"
},

{
    "location": "node_types.html#",
    "page": "Node Types",
    "title": "Node Types",
    "category": "page",
    "text": ""
},

{
    "location": "node_types.html#Node-Types-1",
    "page": "Node Types",
    "title": "Node Types",
    "category": "section",
    "text": "Modules = [PowerDynBase]\nPages   = [\"NodeDynamics/PQAlgebraic.jl\", \"NodeDynamics/PVAlgebraic.jl\", \"NodeDynamics/SlackAlgebraic.jl\", \"NodeDynamics/SwingEquation.jl\"]"
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
    "location": "custom_node_types.html#PowerDynBase.AbstractNodeParameters",
    "page": "Custom Node Types",
    "title": "PowerDynBase.AbstractNodeParameters",
    "category": "type",
    "text": "Abstract super type for all node parameter types.\n\n\n\n\n\n"
},

{
    "location": "custom_node_types.html#Custom-Node-Types-1",
    "page": "Custom Node Types",
    "title": "Custom Node Types",
    "category": "section",
    "text": "To define your own Node Types, use the PowerDynBase.@DynamicNode macro. The new node type will be a subtype of PowerDynBase.AbstractNodeParameters.@DynamicNode\nAbstractNodeParameters"
},

{
    "location": "error_types.html#",
    "page": "Error Types",
    "title": "Error Types",
    "category": "page",
    "text": ""
},

{
    "location": "error_types.html#PowerDynBase.GridDynamicsError",
    "page": "Error Types",
    "title": "PowerDynBase.GridDynamicsError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong during the grid dynamics construction.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynBase.NodeDynamicsError",
    "page": "Error Types",
    "title": "PowerDynBase.NodeDynamicsError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong during the node dynamics construction.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#PowerDynBase.StateError",
    "page": "Error Types",
    "title": "PowerDynBase.StateError",
    "category": "type",
    "text": "Error to be thrown if something goes wrong when creating or modifying states.\n\n\n\n\n\n"
},

{
    "location": "error_types.html#Error-Types-1",
    "page": "Error Types",
    "title": "Error Types",
    "category": "section",
    "text": "Modules = [PowerDynBase]\nPages   = [\"Errors.jl\"]"
},

{
    "location": "internalsBase.html#",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.jl",
    "category": "page",
    "text": ""
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractAlgebraicGridDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractAlgebraicGridDynamics",
    "category": "type",
    "text": "Abstract super type for all grid dynamics represented by DAEs.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractAlgebraicNodeDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractAlgebraicNodeDynamics",
    "category": "type",
    "text": "Abstract super type for all node dynamics represented by DAEs.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractDAEVariable",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractDAEVariable",
    "category": "type",
    "text": "Abstract super type for all Variables for DAE-type node dynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractDEVariable",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractDEVariable",
    "category": "type",
    "text": "Abstract super type for all variables that AbstractNodeDynamics sub types can be called with.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractNetworkFunction",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractNetworkFunction",
    "category": "type",
    "text": "abstract type AbstractNetworkFunction{T<:AbstractNodeDynamics, M<:AbstractMatrix} end\n\nAbstract super type of all functions that define how a differential equation for the whole network / power grid behaves, e.g. the full right-hand-side function of the ODE.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractNodeDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractNodeDynamics",
    "category": "type",
    "text": "Abstract super type for all abstract node dynamics types.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractODEVariable",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractODEVariable",
    "category": "type",
    "text": "Abstract super type for all Variables for ODE-type node dynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractOrdinaryGridDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractOrdinaryGridDynamics",
    "category": "type",
    "text": "Abstract super type for all grid dynamics represented by ODEs.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AbstractOrdinaryNodeDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AbstractOrdinaryNodeDynamics",
    "category": "type",
    "text": "Abstract super type for all node dynamics represented by ODEs.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AlgebraicGridDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AlgebraicGridDynamics",
    "category": "type",
    "text": "TBD\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.AlgebraicNodeDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.AlgebraicNodeDynamics",
    "category": "type",
    "text": "DOCS TBD!\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.BaseState",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.BaseState",
    "category": "type",
    "text": "\n    BaseState(grid, vec)\n\n\nEncode a state vector and the corresponding rhs information.\n\nKeyword Arguments\n\ngrid is a GridDynamics instance that contains the overall system rhs.\nvec is a state vector of the system who\'s length is given by the total       number of internal and voltage variables.\n\nIndexing\n\nIn an instance b of of a BaseState behaves like an Array, i.e. you can access the j-th element of the state vector (and set it to a value ξ) by calling b[j] ( = ξ ).\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.DAEVariable",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.DAEVariable",
    "category": "type",
    "text": "Variables for DAE-type node dynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.ODEVariable",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.ODEVariable",
    "category": "type",
    "text": "Variables for ODE-type node dynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.OrdinaryGridDynamics",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.OrdinaryGridDynamics",
    "category": "type",
    "text": "TBD\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.OrdinaryGridDynamicsWithMass",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.OrdinaryGridDynamicsWithMass",
    "category": "type",
    "text": "TBD\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#Base.view-Tuple{PowerDynBase.AbstractDEVariable,Any}",
    "page": "PowerDynBase.jl",
    "title": "Base.view",
    "category": "method",
    "text": "Extend view from arrays to subtypes of AbstractDEVariable.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.DynamicNode-NTuple{4,Any}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.DynamicNode",
    "category": "method",
    "text": "See DPSABase.@DynamicNode.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase._GridDynamics-Tuple{AbstractArray{#s30,1} where #s30<:OrdinaryNodeDynamics,AbstractArray{T,2} where T}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase._GridDynamics",
    "category": "method",
    "text": "Create for each subtype of DPSABase.AbstractNodeDynamics the corresponding subtype of DPSABase.GridDynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.checkLY-Tuple{AbstractArray{T,2} where T}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.checkLY",
    "category": "method",
    "text": "Check whether the admittance laplacian has no purely nodal admittances, i.e. that the sum of columns and rows equals to zero.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.complexview-Tuple{PowerDynBase.AbstractDEVariable,Any,Any}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.complexview",
    "category": "method",
    "text": "Extend complexview from arrays to subtypes of AbstractDEVariable.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.complexview-Union{Tuple{T}, Tuple{AbstractArray{T,N} where N,Any,Any}} where T",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.complexview",
    "category": "method",
    "text": "Interpret (part of) an array of real values as an array with complex values.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.create_DEVariable-Tuple{Any,Symbol}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.create_DEVariable",
    "category": "method",
    "text": "Basically, this macro generates all the constructors (internal and external) for a subtype of AbstracDEVariable.\n\nIf you really want to understand this macro, uncomment the println(ex) statement at the end and check the resulting generated expression.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.excomparison-Tuple{Any}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.excomparison",
    "category": "method",
    "text": "Create an expresseion where == is applied between all the expressions given as argument here.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.getDEVariableType-Tuple{Type{Val{OrdinaryNodeDynamics}}}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.getDEVariableType",
    "category": "method",
    "text": "Identify each subtype of AbstractNodeDynamics with its corresponding subtype of AbstractDEVariable\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.internal_unitranges-Tuple{AbstractArray{#s30,1} where #s30<:PowerDynBase.AbstractNodeDynamics}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.internal_unitranges",
    "category": "method",
    "text": "Get the unit ranges that indicate where in the array the internal variables for each of the nodes is saved.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.mapfields-Tuple{Any,Any,Vararg{Any,N} where N}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.mapfields",
    "category": "method",
    "text": "function mapfields(f, s, args...)\n\nApplies f to all fields of (the struct) s giving args... as additional arguments.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.nint-Tuple{OrdinaryNodeDynamics}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.nint",
    "category": "method",
    "text": "Get number of internal arguments of the node.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.nodeiterator-Tuple{NetworkRHS,PowerDynBase.AbstractDEVariable,Any}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.nodeiterator",
    "category": "method",
    "text": "nodeiterator(rhs::NetworkRHS, x::AbstractDEVariable, t)\n\nDistribute the values in x over all the nodes that are summarized in rhs.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.rhs2root-Tuple{Function}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.rhs2root",
    "category": "method",
    "text": "A function converting a rhs-type function to a root-type function.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.total_nint-Tuple{AbstractArray{#s30,1} where #s30<:PowerDynBase.AbstractNodeDynamics}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.total_nint",
    "category": "method",
    "text": "Get the total number of internal variables for an array of node dynamics.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.total_nvars-Tuple{AbstractArray{#s30,1} where #s30<:PowerDynBase.AbstractNodeDynamics}",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.total_nvars",
    "category": "method",
    "text": "Get the total number of dynamic variables for an array of node dynamics.\n\nThis is basically the (real) dimension of the system, hence the sum of internal dynamic variables + 2*(number of nodes = number of complex voltages). The 2 is due to the fact that the complex voltages are treated as two real variables.\n\n\n\n\n\n"
},

{
    "location": "internalsBase.html#PowerDynBase.jl-1",
    "page": "PowerDynBase.jl",
    "title": "PowerDynBase.jl",
    "category": "section",
    "text": "Modules = [PowerDynBase]\nPublic = false"
},

{
    "location": "internalsSolve.html#",
    "page": "PowerDynSolve.jl",
    "title": "PowerDynSolve.jl",
    "category": "page",
    "text": ""
},

{
    "location": "internalsSolve.html#PowerDynSolve.jl-1",
    "page": "PowerDynSolve.jl",
    "title": "PowerDynSolve.jl",
    "category": "section",
    "text": "PowerDynSolve.jl is just a helper library for solving a differential equation system created with PowerDynamics.jl. Documentation will come after after stabilization of the api.Modules = [PowerDynSolve]\nPublic = false"
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
    "text": ""
},

{
    "location": "contact.html#Contact-1",
    "page": "Contact",
    "title": "Contact",
    "category": "section",
    "text": "In case of questions simply submit an issue on gitlab. I you don\'t want to contact us publicly, send an email to Tim Kittel (tim.kittel@elena-international.com) or Sabine Auer (sabine.auer@elena-international.com)."
},

]}
