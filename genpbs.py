import numpy
import os

os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib')

global omegam, omegab, omegal, hubble_h, ngrid, boxsize
global execfile, option
global dendir,srcdirbase,logbase,sumdir,inputbase,samdirbase,outputdirbase,zlistfile,z2listfile
global pbsdir
omegam = 0.27
omegab = 0.045
omegal = 0.73
hubble_h = 0.7
ngrid = 306
boxsize = 47.0

option = 2  # To use lgalaxy+semi_rt
execfile = "./ionz_main"

densdir="/research/prace/sph_smooth_cubepm_130315_6_1728_47Mpc_ext2/nc306/"

srcdirbase="srcs/"
samdirbase="sams/"
outputdirbase = "xfrac/"
logbase="logs/"
inputbase = "inputs/"
sumdir="summary/"
pbsdir="pbs/"

zlistfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z3.txt"
z2listfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"


def submit_job(nion):
    pbsfile=pbsdir+"%4.2f.pbs"%(nion)
    nion_list = "nion.list"
    f = open("nion.list","w+")
    print >> f, "1"
    print >> f, nion
    f.close()
    f = open(pbsfile,"w+")
    print >> f, '#!/bin/bash'
    print >> f, '#$ -N LGALAXY_RT'
    print >> f, '#$ -cwd' 
    print >> f, '#$ -pe openmpi 128' 
    print >> f, '#$ -q mps.q@@mps_amd'
    print >> f, '#$ -S /bin/bash'
    
    # source modules environment:
    print >> f, "" 
    print >> f, 'module add sge' 
    print >> f, 'module add gcc/4.8.1' 
    print >> f, 'module add intel-mpi/64/4.1.1/036'
    print >> f, 'module add gsl/gcc/1.15' 

    print >> f, 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib' 
    
    outputdir = outputdirbase+"%4.2f"%(nion)
    os.system("mkdir -p "+outputdir)
    srcdir = srcdirbase+"%4.2f"%(nion)
    os.system("mkdir -p "+srcdir)

    logfile = logbase+"%4.2f"%(nion) # 
    summaryfile = sumdir+"%4.2f"%(nion)+".summary" #

    zf = open(zlistfile,"r")
    z3list = zf.readlines()
    zf.close;

    zf = open(z2listfile,"r")
    z2list = zf.readlines()
    zf.close;

    if len(z3list) != len(z2list):
        print "Error: z2 != z2"

    print >> f, "echo >",logfile
    for i in range(len(z3list)):
        z2 = z2list[i].strip()
        z3 = z3list[i].strip()
        if i == 0:
            prev_z = "-1"
        else:
            prev_z = z3list[i-1].strip()
        denfile = densdir+"/"+z3+"n_all.dat"
        srcfile = srcdir+"/"+z2+".dat"
        print >> f, "echo 'z = "+z3+"'"
        print >> f, 'mpirun -np $NSLOTS numactl -l',execfile,option,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdir,summaryfile, ">>",logfile
        
    print >> f, "echo 'SEQUENCE COMPLETED' >>",logfile
    f.close


submit_job(40000.)
