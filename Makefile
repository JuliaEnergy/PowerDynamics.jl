DOCSPATH=$(HOME)/work/JuliaDPSA/dpsadocs
DOCS_SERVER=dpsadocs
DOCSPATH_SERVER=dpsadocs

.PHONY : docs
deploy-docs: docs push-docs pull-docs-remote

docs:
	julia --color=yes docs/make.jl

push-docs:
	cp -rf docs/build/* $(DOCSPATH)/
	git -C $(DOCSPATH) add $(DOCSPATH)/*
	git -C $(DOCSPATH) commit -m "docs update"
	git -C $(DOCSPATH) push

pull-docs-remote:
	curl https://elena-international.com/update-dpsa-docs.php
#	ssh $(DOCS_SERVER) "cd $(DOCSPATH_SERVER); \
#	       	git stash; \
#		git stash drop; \
#		git pull; \
#		chmod a+rx *" 

		
