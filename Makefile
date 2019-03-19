.PHONY : open-docs docs test coverage clean

ifeq ($(OS),Windows_NT)
    open_cmd=start
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
    	open_cmd=xdg-open
    else ifeq ($(UNAME_S),Darwin)
    	open_cmd=open
    endif
endif

open-docs: docs
	$(open_cmd) docs/build/index.html

docs:
	julia -e "using Pkg; Pkg.activate(\".\"); include(\"docs/make.jl\")"

test:
	julia --check-bounds=yes --color=yes -e "using Pkg; Pkg.activate(\".\"); Pkg.test()"

coverage:
	julia --check-bounds=yes --color=yes --inline=no -e "using Pkg; Pkg.activate(\".\"); Pkg.test(coverage=true)"

clean:
	find src -name "*.jl.*.cov" -type f -delete
	find test -name "*.jl.*.cov" -type f -delete
	rm -rf docs/build
