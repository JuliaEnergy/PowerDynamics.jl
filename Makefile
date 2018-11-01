.PHONY : docs test coverage clean

docs:
	julia -e "using Pkg; Pkg.activate(\".\"); include(\"docs/make.jl\")"

test:
	julia --check-bounds=yes --color=yes -e "using Pkg; Pkg.activate(\".\"); Pkg.test()"

coverage:
	julia --check-bounds=yes --color=yes --inline=no -e "using Pkg; Pkg.activate(\".\"); Pkg.test(coverage=true)"

clean:
	find src -name "*.jl.*.cov" -type f -delete
	find test -name "*.jl.*.cov" -type f -delete
