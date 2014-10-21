import numpy
import os

os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib')

#lgalaxy things
global firstfile,lastfile,SimulationDir, maxmemsize, lgal_template
global lgal_exec
lgal_exec = "/mnt/lustre/scratch/cs390/codes/47Mpc/L-Galaxies_development/L-Galaxies"
firstfile = 0
lastfile = 127
SimulationDir = "/mnt/lustre/scratch/cs390/47Mpc/"
maxmemsize = 4000
template = "/mnt/lustre/scratch/cs390/47Mpc/couple/model_001/input_47mpc_template.par"
zout_list = "zout_lgal"

lgal_template = open(template,"r").read()


#create sources things
global gensource

gensource_exec = "/mnt/lustre/scratch/cs390/codes/Sources_generate/gensource"
#ionz things
global omegam, omegab, omegal, hubble_h, ngrid, boxsize
global ionz_execfile, option
global dendir,srcdirbase,logbase,sumdir,inputbase,samdirbase,outputdirbase,zlistfile,z2listfile
global pbsdir

omegam = 0.27
omegab = 0.045
omegal = 0.73
hubble_h = 0.7
ngrid = 306
boxsize = 47.0

option = 2  
ionz_execfile = "/mnt/lustre/scratch/cs390/codes/ionz_codes/ionz_main"

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


os.system("mkdir -p "+logbase)
os.system("mkdir -p "+inputbase)
os.system("mkdir -p "+sumdir)
os.system("mkdir -p "+pbsdir)
os.system("mkdir -p "+srcdirbase)
os.system("mkdir -p "+samdirbase)
os.system("mkdir -p "+outputdirbase)

def submit_job(nion):
    # prepare PBS file
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
    
    # prepare folder/files
    outputdir = outputdirbase+"%4.2f"%(nion)
    os.system("mkdir -p "+outputdir)
    srcdir = srcdirbase+"%4.2f"%(nion)
    os.system("mkdir -p "+srcdir)
    
    samdir = samdirbase+"%4.2f"%(nion)
    os.system("mkdir -p "+samdir)
    
    zout_file = samdirbase+"z_list"
    logfile = logbase+"%4.2f"%(nion)+".ionz" # 
    lgal_log = logbase+"%4.2f"%(nion)+".lgal" # 
    convert_log = logbase+"%4.2f"%(nion)+".convert" # 
    summaryfile = sumdir+"%4.2f"%(nion)+".sum" #
    lgal_input = inputbase+"lgal_%4.2f"%(nion)
    # prepare input file for L-galaxy
    lgalinp = open(lgal_input,"w")
    print >> lgalinp, "FirstFile",firstfile
    print >> lgalinp, "LastFile", lastfile
    print >> lgalinp, "OutputDir", samdir
    print >> lgalinp, "SimulationDir", SimulationDir
    print >> lgalinp, "ReionizationOn  5"
    print >> lgalinp, "Reionization_z0  8.0"
    print >> lgalinp, "Reionization_zr  7.0"
    print >> lgalinp, "XfracDir", outputdir
    print >> lgalinp, "FileWithOutputRedshifts", zout_file 
    print >> lgalinp, "MaxMemSize" ,maxmemsize
    print >> lgalinp, lgal_template
    lgalinp.close()

    zf = open(zlistfile,"r")
    z3list = zf.readlines()
    zf.close;

    zf = open(z2listfile,"r")
    z2list = zf.readlines()
    zf.close;

    if len(z3list) != len(z2list):
        print "Error: z2 != z2"
        exit()
    print >> f, "rm -f "+lgal_log
    print >> f, "rm -f "+convert_log
    print >> f, "rm -f "+logfile

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
        print >> f, "echo '",z2,"' > ",zout_file
        print >> f, "# run lgalaxy"
        print >> f, 'mpirun -np $NSLOTS numactl -l',lgal_exec,lgal_input,">>",lgal_log
        print >> f, '# run gensourc for current snapshot'
        print >> f, gensource_exec,i,samdir+"/SA_z",srcdir,z2listfile,">>", convert_log 
        print >> f, '# run ionz'
        print >> f, 'mpirun -np $NSLOTS numactl -l',ionz_execfile,option,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdir,summaryfile, ">>",logfile
        
    print >> f, "echo '#SEQUENCE COMPLETED' >>",summaryfile
    f.close


submit_job(40000.)
