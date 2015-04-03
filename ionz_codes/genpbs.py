import numpy
pbsfile="runall.pbs"
f = open(pbsfile,"w+")
print >> f, '#!/bin/bash'
print >> f, '#$ -N custom'
print >> f, '#$ -cwd' 
print >> f, '#$ -pe openmpi 128' 
print >> f, '#$ -q mps.q'
print >> f, '#$ -S /bin/bash'

# source modules environment:
print >> f, "" 
print >> f, 'module add sge' 
print >> f, 'module add gcc/4.8.1' 
print >> f, 'module add intel-mpi/64/4.1.1/036'
print >> f, 'module add gsl/gcc/1.15' 

print >> f, 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib' 

#parameters
execfile = "./ionz_main"
nion_list = "nion.list"
omegam = 0.27
omegab = 0.044
omegal = 0.73
hubble_h = 0.7
ngrid = 306
boxsize = 47.0

densdir="/research/prace/sph_smooth_cubepm_130315_6_1728_47Mpc_ext2/nc306/"
srcdir="/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/sources/"
outputdir = "./"
zlistfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z3.txt"
z2listfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"


zf = open(zlistfile,"r")
z3list = zf.readlines()
zf.close;

zf = open(z2listfile,"r")
z2list = zf.readlines()
zf.close;

if len(z3list) != len(z2list):
    print "Error: z2 != z2"

print >> f, "echo > log"
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
    print >> f, 'mpirun -np $NSLOTS',execfile,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdir,">> log"

#./ionz_main nion.list 0.27 0.044 0.73 0.7 306 47.0 /research/prace/sph_smooth_cubepm_130315_6_1728_47Mpc_ext2/nc306/6.000n_all.dat /mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/sources/6.00.dat 6.000 6.056 outputfolder
f.close
