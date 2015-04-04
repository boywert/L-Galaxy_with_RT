all: genpbs ionz

genpbs: genpbs.py config.py
	python genpbs.py
ionz:
	cd ionz_codes
	make
	cd ..
clean:
	rm -f *~ log/* pbs/* 
