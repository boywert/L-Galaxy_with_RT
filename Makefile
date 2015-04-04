reionization_code="ionz_codes/"
gen_source_code="sources_generate/"
sam_code="L-Galaxies_development/"
lgal_struct="L-Galaxies_development/awk/L-Galaxies.h"

all: genpbs reion sam gen_source

genpbs: genpbs.py config.py
	python genpbs.py
reion:
	cd $(reionization_code) && $(MAKE) clean && $(MAKE)
	cd ..
sam:
	cd $(sam_code) && $(MAKE) clean && $(MAKE) && $(MAKE) metadata
	cd ..
gen_source: sam
	cp $(lgal_struct) $(gen_source_code)
	cd $(gen_source_code) && $(MAKE) clean && $(MAKE)
	cd ..
clean:
	rm -f *~ log/* pbs/* 
