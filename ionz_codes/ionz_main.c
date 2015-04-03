/**
 * @file   ionz_main.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 21:01:59 2014
 * 
 * @brief  Main program
 * 
 */

#include "ion.h"
struct_const constvars = {3.14159265359,1024,0.1,2.0};
/** 
 * Main program
 * 
 * @param argc 
 * @param argv 
 * 
 * @return 
 */
int main(int argc, char **argv) {
  int ii, jj, jk;
  float r_min,r_max;
  int Nnion,N1,N2,N3;
  float *nion; 
  double robar,robarhalo;
  float *Radii_list;
  int n_radii;
  int *NjobsperTask;
  int *JobsTask;
  double t_start, t_stop;
  float *buffer, *buffer_final;
  float vomegam,vomegab,vomegalam;
  fftw_real ***nh, ***ngamma, ****nxion;
  fftw_real ****xfrac;
#ifdef CHUNKTRANSFER
  int mpi_buffer=1000000;
  int cur_len;
#endif
  char inputfile[2000];
  char densfilename[2000], sourcefilename[2000];
  char z_prev[1000],z_out[1000];
  char outputdir[2000];


  int use_prev_xfrac = 0;


#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mympi.ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &mympi.NTask);
#else
  mympi.NTask = 1;
  mympi.ThisTask = 0;
#endif //PARALLEL

  if(mympi.ThisTask == 0) {
    if(argc == 2) {
      sprintf(inputfile, "%s", argv[1]);
      printf("Read all parameters from input file: %s\n",inputfile);
      printf("I still don't have time to do this option....Boyd\n");
      exit(1);
    }
    else if(argc > 2 ) {
      printf("Arguments:\n");
      printf("option:\t%s\n",argv[1]);
      sscanf(argv[1],"%d",&input_param.option);
      printf("nion_file:\t%s\n",argv[2]);
      read_nion(argv[2]);
      printf("Omega_matter:\t%s\n",argv[3]);
      sscanf(argv[3],"%f",&input_param.omegam);
      printf("Omega_baryon:\t%s\n",argv[4]);
      sscanf(argv[4],"%f",&input_param.omegab);
      printf("Omega_Lambda:\t%s\n",argv[5]);
      sscanf(argv[5],"%f",&input_param.omegalam);
      printf("Hubble_h:\t%s\n",argv[6]);
      sscanf(argv[6],"%f",&input_param.Hubble_h);
      printf("N_grid:\t%s\n",argv[7]);
      sscanf(argv[7],"%d",&input_param.N1);
      input_param.N2 = input_param.N1;
      input_param.N3 = input_param.N1;
      printf("Boxsize:\t%s\n",argv[8]);
      sscanf(argv[8],"%f",&input_param.boxsize);
      input_param.gridsize = input_param.boxsize/(float)input_param.N1;
      printf("densityfile:\t%s\n",argv[9]);
      sprintf(input_param.densityfile,"%s",argv[9]);
      printf("sourcesfile:\t%s\n",argv[10]);
      sprintf(input_param.sourcesfile,"%s",argv[10]);
      printf("Current redshift:\t%s\n",argv[11]);
      sprintf(input_param.cur_z,"%s",argv[11]);
      printf("Previous redshift:\t%s\n",argv[12]);
      sprintf(input_param.prev_z,"%s",argv[12]);
      printf("Output folder:\t%s\n",argv[13]);
      sprintf(input_param.outputdir,"%s",argv[13]);
      printf("Summary file:\t%s\n",argv[14]);
      sprintf(input_param.summary_file,"%s",argv[14]);
                            
    }
    else {
      printf("Usage[1]: ./exec inputfile\n");
      printf("Usage[2]: ./exec option nion_file omegam omegab omegalam hubble_h n_grid boxsize densityfile sourcefile curr_z prev_z outputfolder (very useful for submitting batch MPI tasks)\n");
      exit(1);
    }
    if(mympi.ThisTask == 0) {
      printf("Start semi-numerical reionization process\n");
    }  
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);    
  MPI_Bcast(&input_param.option, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.Nnion, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.nion[0], 100, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.a_expansion, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.z, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.Hubble_h, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.omegam, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.omegalam, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.omegab, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.N1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.N2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.N3, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.boxsize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.gridsize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.densityfile[0],2000, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.sourcesfile[0],2000, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.cur_z[0],100, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.prev_z[0],100, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.outputdir[0],2000, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_param.summary_file[0],2000, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  Nnion = input_param.Nnion;
  nion=(float*)calloc(Nnion,sizeof(float));
  for(ii=0;ii<Nnion;ii++) {
    nion[ii] = input_param.nion[ii];
  }

  N1 = input_param.N1;
  N2 = input_param.N2;
  N3 = input_param.N3;
  vomegam = input_param.omegam;
  vomegalam = input_param.omegalam;
  vomegab = input_param.omegab;
  if(input_param.option == 2)
    use_prev_xfrac = 1;
  sprintf(densfilename,"%s",input_param.densityfile);
  sprintf(sourcefilename,"%s",input_param.sourcesfile);
  sprintf(z_out,"%s",input_param.cur_z);
  sprintf(z_prev,"%s",input_param.prev_z);
  sprintf(outputdir,"%s",input_param.outputdir);
  if(mympi.ThisTask == 0) {
    printf("Using Cosmological parameters:\t");
    printf("Omega_m:%f\t",vomegam);
    printf("Omega_b:%f\t",vomegab);
    printf("Omega_lambda:%f\n",vomegalam);
    printf("Grid:\t%dx%dx%d\n",N1,N2,N3);
  }

  /* Allocating memory to different arrays */
  // Setting_Up_Memory_For_ionz(Nnion, N1, N2, N3
  nh = allocate_fftw_real_3d(N1,N2,N3+2);
  ngamma = allocate_fftw_real_3d(N1,N2,N3+2);
  nxion=(fftw_real****)malloc(sizeof(fftw_real***)*Nnion);
  if(use_prev_xfrac == 1)
    xfrac=(fftw_real****)malloc(sizeof(fftw_real***)*Nnion);
  for(jk=0;jk<Nnion;++jk) {
    nxion[jk] = allocate_fftw_real_3d(N1,N2,N3+2);
    if(use_prev_xfrac == 1)
      xfrac[jk] = allocate_fftw_real_3d(N1,N2,N3+2);
  }
  t_start =Get_Current_time();
  
  /* Allocate buffer to store 3D array */
  buffer = malloc(sizeof(float)*N1*N2*N3);
  /* Use Task:0 to read density */
  if(mympi.ThisTask == 0) {
    read_density(densfilename,buffer,&robar,N1,N2,N3,vomegam,vomegab);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  chunk_float_mpi_bcast(buffer,N1*N2*N3,mpi_buffer,0);
#else 
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif 
  MPI_Barrier(MPI_COMM_WORLD);
#endif // PARALLEL
  
  unpack_3d_array_mpi_transfer(buffer,nh,N1,N2,N3);
  
  
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(mympi.ThisTask == 0)  {
    read_sources(sourcefilename,buffer,&robarhalo,N1,N2,N3);  
  } 

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robarhalo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  chunk_float_mpi_bcast(buffer,N1*N2*N3,mpi_buffer,0);
#else 
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif 
#endif
  unpack_3d_array_mpi_transfer(buffer,ngamma,N1,N2,N3);
  free(buffer);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(mympi.ThisTask == 0) 
    printf("roHIbar = %lf, roUV = %lf\n",robar,robarhalo);


  /* Read xfrac from previous snapshot */
  if(use_prev_xfrac == 1) {
    buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
    if(mympi.ThisTask == 0)
      read_xfrac(outputdir, z_prev, buffer, nion, Nnion, N1, N2, N3);
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
    chunk_float_mpi_bcast(buffer,Nnion*N1*N2*N3,mpi_buffer,0);
#else 
    MPI_Bcast(buffer, Nnion*N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif 
#endif

    unpack_4d_array_mpi_transfer(buffer, xfrac, Nnion, N1, N2, N3);

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  t_stop = Get_Current_time();

#ifdef PARALLEL
  /* Sanity MPI check */
  if(mympi.ThisTask == 1)
    printf("N1=%d N2=%d N3=%d\n",N1,N2,N3);
#endif
  if(mympi.ThisTask == 0)
    printf("reading in data %lf s\n",t_stop-t_start);

  //calculating max and min radius for smoothing in units of grid size
  r_min=1.;
  r_max=pow((1.*N1*N2*N3),(1./3.))/2.;

  Radii_list = malloc(sizeof(float)*constvars.max_Nradii); 
  NjobsperTask = malloc(sizeof(float)*mympi.NTask);
  n_radii = make_radii_list(Radii_list,r_min,r_max,constvars.dr_inc,constvars.max_dr);
  for(jj=0;jj<mympi.NTask;jj++) {
    NjobsperTask[jj] = n_radii/mympi.NTask;
    if(jj < n_radii%mympi.NTask)
      NjobsperTask[jj]++;
    if(jj == mympi.ThisTask) {
      JobsTask = malloc(sizeof(int)*NjobsperTask[jj]);
      for(ii=0;ii<NjobsperTask[jj];ii++) 
	JobsTask[ii] = ii*mympi.NTask+mympi.ThisTask;
    }
  }


#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif  
  t_start = Get_Current_time();
  // Do subgrid seminumerical simulation
  if(mympi.ThisTask == 0)
    printf("Start subgrid semi-numerical reionization\n");
  if(use_prev_xfrac == 1)
    subgrid_reionization_with_xfrac(nh, ngamma, xfrac, nxion, robar, nion, Nnion, N1, N2, N3 );  
  else
    subgrid_reionization(nh, ngamma, nxion, robar, nion, Nnion, N1, N2, N3 );
 
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_stop = Get_Current_time();
  if(mympi.ThisTask == 0)
    printf("Finish subgrid semi-numerical reionization: %lf s\n",t_stop-t_start);
  if(mympi.ThisTask == 0)
    printf("Start semi-numerical reionization process\n");

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_start = Get_Current_time();
  for(ii=0;ii<NjobsperTask[mympi.ThisTask];ii++) { 
    if(use_prev_xfrac == 1)
      reionization_with_xfrac(Radii_list[JobsTask[ii]], nh, ngamma, xfrac, nxion, nion, Nnion, N1, N2, N3);
    else
      reionization(Radii_list[JobsTask[ii]], nh, ngamma, nxion, nion, Nnion, N1, N2, N3);
  }
  free_fftw_real_3d(ngamma,N1,N2,N3+2);
  if(use_prev_xfrac == 1) {
    for(ii=0;ii<Nnion;ii++)
      free_fftw_real_3d(xfrac[ii],N1,N2,N3+2);
    free(xfrac);
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
  t_stop = Get_Current_time();
  if(mympi.ThisTask == 0)
    printf("Finish reionizing process %lf s\n",t_stop-t_start);

  t_start = Get_Current_time();
  pack_4d_array_mpi_transfer(nxion,buffer,Nnion, N1, N2, N3);
  for(ii=0;ii<Nnion;ii++)
    free_fftw_real_3d(nxion[ii],N1,N2,N3+2);
  free(nxion);
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  t_stop = Get_Current_time();
  if(mympi.ThisTask == 0) {
    printf("Finish packing data %lf s\n",t_stop-t_start);
#ifdef PARALLEL
    buffer_final = malloc(sizeof(float)*Nnion*N1*N2*N3);
#endif
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  t_start = Get_Current_time();
  ii = 0;
  while (ii*mpi_buffer < Nnion*N1*N2*N3) {
    cur_len = min(mpi_buffer, Nnion*N1*N2*N3-ii*mpi_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&buffer[ii*mpi_buffer],&buffer_final[ii*mpi_buffer],cur_len,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
    ii++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = Get_Current_time();
  if(mympi.ThisTask == 0)
    printf("Finish finding max:split %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#else // ~CHUNKTRANSFER
  t_start = Get_Current_time();
  MPI_Reduce(buffer, buffer_final, Nnion*N1*N2*N3, MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = Get_Current_time();
  if(mympi.ThisTask == 0)
    printf("Finish finding max:whole %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#endif // CHUNKTRANSFER
#endif // PARALLEL  
  if(mympi.ThisTask == 0) {
#ifdef PARALLEL
    write_xfrac(outputdir, z_out, buffer_final, nh, robar,nion, Nnion, N1, N2, N3);
#else
    write_xfrac(outputdir, z_out, buffer, nh, nion, robar, Nnion, N1, N2, N3);
#endif
  }
  MPI_Finalize();
  return 0;
} /* END OF MAIN */

