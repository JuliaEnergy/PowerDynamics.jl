module PowerDynamics

using Reexport: Reexport, @reexport
@reexport using NetworkDynamics
using NetworkDynamics: SymbolicView, getcomp

using SciMLBase: SciMLBase
using NonlinearSolve: NonlinearSolve
using ForwardDiff: ForwardDiff
using LinearAlgebra: LinearAlgebra
using Graphs: nv, ne
using ArgCheck: @argcheck
using Setfield: @set, @set!
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using MacroTools: postwalk, @capture

using ModelingToolkit: ModelingToolkit, @connector, @named,
                       System, getname, unknowns, get_name, t_nounits as t, Equation,
                       mtkcompile
# needed for @mtkmodel
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect
using Symbolics: Symbolics
using SciMLBase: SciMLBase

# imports to load sparsity extensions
using SparseConnectivityTracer: SparseConnectivityTracer
using SparseMatrixColorings: SparseMatrixColorings

export Terminal, BusBar, LineEnd
export MTKBus, MTKLine, CompositeInjector, Ibase, Zbase, Ybase
include("modeling_tools.jl")

export isinjectormodel, isbusmodel, isbranchmodel, islinemodel
include("interfaces.jl")

export @attach_metadata!, set_voltage!, set_current!, refine_timeseries
include("utils.jl")
include("Library/Library.jl")
using .Library

export compile_line, compile_bus
export Line, Bus
export simplify_mtkline, simplify_mtkbus
include("network_components.jl")

using DataFrames: DataFrame
using OrderedCollections: OrderedDict
export pfSlack, pfPV, pfPQ
export solve_powerflow, initialize_from_pf!, initialize_from_pf, show_powerflow
export powerflow_model, ispfmodel
export has_pfmodel, get_pfmodel, set_pfmodel!, delete_pfmodel!
include("powerflow.jl")

export PFInitConstraint, @pfinitconstraint, PFInitFormula, @pfinitformula, copy_pf_parameters
export add_pfinitconstraint!, add_pfinitformula!
export set_pfinitconstraint!, set_pfinitformula!
export has_pfinitconstraint, has_pfinitformula
export get_pfinitconstraints, get_pfinitformulas
export delete_pfinitconstraints!, delete_pfinitformulas!
include("pfinitconstraint.jl")

end
