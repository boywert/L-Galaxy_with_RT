This_dir=$(CURDIR)
gensrc_dir=“sources_generate”

all: sourcegen

sourcegen:
	cd $(gensrc_dir)
	make clean all
	cd $(This_dir)
