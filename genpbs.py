import numpy
import os
import re
import config

this_dir = os.getcwd()
os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib')

#lgalaxy things
global firstfile,lastfile,SimulationDir, maxmemsize, lgal_template
global lgal_exec
global lgal_folder

lgal_folder = config.lgal_folder
lgal_exec = lgal_folder+"/L-Galaxies"

firstfile = config.firstfile
lastfile = config.lastfile
SimulationDir = config.SimulationDir
maxmemsize = config.maxmemsize
template = config.template
reionization_option = config.reionization_option

lgal_template = open(template,"r").read()
lgal_struct = config.lgal_struct

#create sources things
global gensource
gensource_dir = config.gensource_dir
gensource_exec = gensource_dir+"/gen_source"
#ionz things
global omegam, omegab, omegal, hubble_h, ngrid, boxsize
global ionz_execfile, option
global dendir,srcdirbase,logbase,sumdir,inputbase,samdirbase,outputdirbase,zlistfile,z2listfile
global pbsdir

omegam = config.omegam
omegab = config.omegab
omegal = config.omegal
hubble_h = config.hubble_h
ngrid = config.ngrid
boxsize = config.boxsize

option = config.option
ionz_execfile = config.ionz_execfile

densdir=config.densdir

srcdirbase=config.srcdirbase
samdirbase=config.samdirbase
outputdirbase = config.outputdirbase
logbase=config.logbase
inputbase = config.inputbase
sumdir=config.sumdir
pbsdir=config.pbsdir

alistfile = config.alistfile
zlistfile=config.zlistfile
z2listfile=config.z2listfile


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
    summaryfile = sumdir+"%4.2f"%(nion)+".sum" 
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
    print >> f, '#$ -q mps.q@@compute_amd_c6145_mps'
    print >> f, '#$ -S /bin/bash'
    print >> f, "#$ -j y"
    print >> f, "#$ -o logs/%4.2f.log"%(nion)
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
    os.system("cp "+lgal_folder+"/python/LGalaxyStruct.py "+samdir)
    zout_file = this_dir+"/"+samdirbase+"z_list"
    logfile = logbase+"%4.2f"%(nion)+".ionz" # 
    lgal_log = logbase+"%4.2f"%(nion)+".lgal" # 
    convert_log = logbase+"%4.2f"%(nion)+".convert" # 
    lgal_input = inputbase+"lgal_%4.2f"%(nion)
    # prepare input file for L-galaxy
    lgalinp = open(lgal_input,"w")
    print >> lgalinp, "FirstFile",firstfile
    print >> lgalinp, "LastFile", lastfile
    print >> lgalinp, "FileWithZList", alistfile
    print >> lgalinp, "OutputDir", samdir
    print >> lgalinp, "SimulationDir", SimulationDir
    print >> lgalinp, "ReionizationOn", reionization_option
    print >> lgalinp, "Reionization_z0  8.0"
    print >> lgalinp, "Reionization_zr  7.0"
    print >> lgalinp, "XfracDir", outputdir
    print >> lgalinp, "XfracNGrids", ngrid
    print >> lgalinp, "FileWithOutputRedshifts", zout_file
    print >> lgalinp, "FileWithZList_OriginalCosm", alistfile
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
        if(i < len(z3list)-1):
            z2_next = z2list[i+1].strip()
        z3 = z3list[i].strip()
        if i == 0:
            prev_z = "-1"
        else:
            prev_z = z3list[i-1].strip()
        denfile = densdir+"/"+z3+"ntot_all.dat"
        srcfile = srcdir+"/"+z2+".dat"
        print >> f, "echo 'z = "+z3+"'"
        print >> f, "echo '"+z2+"' > ",zout_file
        if(i < len(z3list)-1):
            print >> f, "echo '"+z2_next+"' >> ",zout_file
        else:
            print >> f, "echo 30.00 >> ",zout_file
        print >> f, "cat ",zout_file
        print >> f, "# run lgalaxy"
        print >> f, 'echo','mpirun -np $NSLOTS numactl -l',lgal_exec,lgal_input
        print >> f, 'time mpirun -np $NSLOTS numactl -l',lgal_exec,lgal_input,">>",lgal_log
        print >> f, '# run gensourc for current snapshot'
        print >> f, 'echo',gensource_exec,i,samdir+"/SA_z",srcdir,z2listfile,ngrid
        print >> f, 'time',gensource_exec,i,samdir+"/SA_z",srcdir,z2listfile,ngrid,">>", convert_log 
        print >> f, '# run ionz'
        print >> f, "echo ", 'mpirun -np $NSLOTS numactl -l',ionz_execfile,option,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdirbase,summaryfile
        print >> f, 'time mpirun -np $NSLOTS numactl -l',ionz_execfile,option,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdirbase,summaryfile, ">>",logfile
        
    print >> f, "echo '#SEQUENCE COMPLETED' >>",summaryfile
    f.close


submit_job(480.)
