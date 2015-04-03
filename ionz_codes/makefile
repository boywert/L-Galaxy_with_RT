export TACC_FFTW2_LIB=/home/c/cs/cs390/local/fftw-2.1.5/install/lib/
export TACC_FFTW2_INC=/home/c/cs/cs390/local/fftw-2.1.5/install/include/
LINKLIB= -L${TACC_FFTW2_LIB} -lsrfftw -lsfftw 
INCLUDE= -I${TACC_FFTW2_INC}
CFLAGS=-g -Wall -std=c99
CC=mpicc
CFLAGS+= -DPARALLEL -DCHUNKTRANSFER -DUSE_FORTRAN_SPEEDUP_ARRAY
POSTFLAGS= -lm

FC=mpif90


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib

all: ionz_main

arrayoperations.o: arrayoperations.f90
	$(FC) -c arrayoperations.f90

read_param.o: read_param.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) read_param.c $(POSTFLAGS)

ionz_misc.o: ionz_misc.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_misc.c $(POSTFLAGS)


ionz_main.o: 	ionz_main.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_main.c $(POSTFLAGS)

ionz_io.o:  ionz_io.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_io.c $(POSTFLAGS)

ionz_funcs.o:	ionz_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_funcs.c $(POSTFLAGS) 


ionz_main: ionz_main.o ionz_misc.o ionz_io.o ionz_funcs.o read_param.o arrayoperations.o
	$(CC) $(CFLAGS) $(INCLUDE) $(LINKLIB) -o ionz_main ionz_main.o ionz_io.o ionz_misc.o ionz_funcs.o read_param.o arrayoperations.o $(POSTFLAGS) 

clean:
	rm -rf *.o
	rm -rf *~








