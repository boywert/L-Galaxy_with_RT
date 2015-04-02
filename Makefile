all: genbps

genpbs: genpbs.py config.py
	python genpbs.py

clean:
	rm -f *~ log/* pbs/* 