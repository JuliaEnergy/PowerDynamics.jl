.PHONY : docs

docs:
	julia -e "using Pkg; Pkg.activate(\".\"); include(\"docs/make.jl\")"

clean:
	find src -name "*.jl.*.cov" -type f -delete
	find test -name "*.jl.*.cov" -type f -delete
