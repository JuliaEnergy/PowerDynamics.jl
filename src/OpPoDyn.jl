module OpPoDyn

using NetworkDynamics: NetworkDynamics, VertexModel, EdgeModel,
                       has_default, set_default!, get_default,
                       has_init, get_init, has_guess, get_guess,
                       has_graphelement, set_graphelement!, get_graphelement, Network,
                       has_metadata, get_metadata, set_metadata!, NWState, uflat, pflat,
                       initialize_componentwise, initialize_componentwise!, interface_values,
                       find_fixpoint, vidxs, extract_nw, SymbolicView,
                       VIndex, EIndex, InitFormula, InitConstraint,
                       sym, psym
using SciMLBase: SciMLBase, solve
using NonlinearSolve: NonlinearSolve, NonlinearProblem
using ForwardDiff: ForwardDiff
using LinearAlgebra: LinearAlgebra, Diagonal, diag, pinv, eigvals
using Graphs: SimpleGraph, add_edge!, nv, ne
using ArgCheck: @argcheck
using Setfield: @set, @set!
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using MacroTools: postwalk, @capture

using ModelingToolkit: ModelingToolkit, @connector, @mtkmodel, @variables, @named,
                       ODESystem, connect, getname, unknowns, get_name, get_iv, get_systems,
                       get_gui_metadata, t_nounits as t, Equation,
                       defaults, parameters, iscomplete, rename, simplify, unwrap,
                       get_eqs, get_observed, get_ctrls, get_defaults, get_schedule,
                       get_connector_type, get_preface, get_initializesystem,
                       get_continuous_events, get_discrete_events, get_parameter_dependencies,
                       get_discrete_subsystems, get_solved_unknowns, get_tspan, get_guesses,
                       structural_simplify
using ModelingToolkit: @unpack, Num, System # needed for @mtkmodel?
using Symbolics: Symbolics, Symbolic, iscall, fixpoint_sub
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve

export Terminal, BusBar, LineEnd
export MTKBus, MTKLine, CompositeInjector, Ibase, Zbase, Ybase
include("modeling_tools.jl")

export isinjectormodel, isbusmodel, isbranchmodel, islinemodel
include("interfaces.jl")

export @attach_metadata!, set_voltage!, set_current!
include("utils.jl")
include("Library/Library.jl")
using .Library

export Line, Bus
export simplify_mtkline, simplify_mtkbus
include("network_components.jl")

using DataFrames: DataFrame
using OrderedCollections: OrderedDict
export pfSlack, pfPV, pfPQ
export solve_powerflow, initialize_from_pf!, initialize_from_pf, show_powerflow
export powerflow_model, ispfmodel
include("powerflow.jl")

export PFInitConstraint, @pfinitconstraint, PFInitFormula, @pfinitformula, copy_pf_parameters
export add_pfinitconstraint!, add_pfinitformula!
export set_pfinitconstraint!, set_pfinitformula!
export has_pfinitconstraint, has_pfinitformula
export get_pfinitconstraints, get_pfinitformulas
export delete_pfinitconstraints!, delete_pfinitformulas!
include("pfinitconstraint.jl")

include("pin_parameters.jl")

end
