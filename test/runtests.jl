using PowerDynamics
using ParallelTestRunner
using Test

# `runtests` runs every entry of `testsuite` in its own module on a worker process. The files are
# listed explicitly instead of autodiscovered, so this file stays the single place which decides
# what runs.

TESTDIR = pkgdir(PowerDynamics, "test")
testsuite = Dict{String,Expr}()

# init_worker_code -> workers Main, init_code -> every test module
init_worker_code = quote
    using PowerDynamics
    PowerDynamics.load_pdtesting()
end

"""
    includetest!(file)

Add a test file:
- plain repl: `include(file)` in Main
- ] test/include("runtests.jl"): add to parallel test pool and run on worker
"""
function includetest!(file)
    path = joinpath(TESTDIR, file)
    if any(f -> f.func === :include || f.func === :_include, stacktrace())
        testsuite[file] = :(Main.PowerDynamicsTesting.quiet_include(@__MODULE__, $path))
    else
        Core.eval(Main, init_worker_code) # no worker to set up here
        @testset "$file" begin include(path) end
    end
end

includetest!("quality_test.jl");
includetest!("Library_test.jl");
includetest!("saturation_test.jl");
includetest!("modeling_tools_test.jl");
includetest!("initialization_test.jl");

# OpenIPSL Ref Test
includetest!("OpenIPSL_test/PSSE_GENCLS_test.jl");
includetest!("OpenIPSL_test/PSSE_GENROU_test.jl");
includetest!("OpenIPSL_test/PSSE_GENROE_test.jl");
includetest!("OpenIPSL_test/PSSE_GENSAL_test.jl");
includetest!("OpenIPSL_test/PSSE_GENSAE_test.jl");
includetest!("OpenIPSL_test/PSSE_IEEET1_test.jl");
includetest!("OpenIPSL_test/PSSE_SCRX_test.jl");
includetest!("OpenIPSL_test/PSSE_ESST4B_test.jl");
includetest!("OpenIPSL_test/PSSE_EXST1_test.jl");
includetest!("OpenIPSL_test/PSSE_ESST1A_test.jl");
includetest!("OpenIPSL_test/PSSE_IEEEG1_test.jl");
includetest!("OpenIPSL_test/PSSE_GGOV1_test.jl");
includetest!("OpenIPSL_test/PSSE_HYGOV_test.jl");
includetest!("OpenIPSL_test/PSSE_IEEEST_test.jl");

includetest!("validation/ieee39_RMSPowerSims.jl/ieee39_validation.jl");

for dir in ("tutorials", "examples")
    for file in readdir(joinpath(TESTDIR, "..", "docs", dir))
        endswith(file, ".jl") || continue
        includetest!("../docs/$dir/" * file);
    end
end

# empty when the lines above were sent to a REPL one by one, those tests already ran
isempty(testsuite) || runtests(PowerDynamics, ARGS; testsuite, init_worker_code)
