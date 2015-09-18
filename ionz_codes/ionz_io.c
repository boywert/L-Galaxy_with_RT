/**
 * @file   ionz_io.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 21:44:53 2014
 * 
 * @brief  IO operations
 * 
 * 
 */

#include "ion.h"

/** 
 * Convert 3D-fftw_real array to 1D-float array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 3D-fftw_real array input
 * @param output 1D-float array output
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void pack_3d_array_mpi_transfer(fftw_real ***input, float *output, int N1, int N2, int N3) { 
  int ii,jj,kk;
  for(kk=0;kk<N3;kk++)
    for(jj=0;jj<N2;jj++)
      for(ii=0;ii<N1;ii++)
	output[kk*N2*N1 + jj*N1 + ii] = input[ii][jj][kk];
}

/** 
 * Convert 1D-float array to 3D-fftw_real array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 1D-float array input
 * @param output 3D-fftw_real array output
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void unpack_3d_array_mpi_transfer(float *input, fftw_real ***output, int N1, int N2, int N3) {
  int ii,jj,kk;
  for(kk=0;kk<N3;kk++)
    for(jj=0;jj<N2;jj++)
      for(ii=0;ii<N1;ii++) { 
	// if(mympi.ThisTask == 0) printf("%d %d %d %f\n",ii,jj,kk,output[ii][jj][kk]);
	// if(mympi.ThisTask == 0) printf("%d %f\n",kk*N2*N1 + jj*N1 + ii,input[kk*N2*N1 + jj*N1 + ii]);
	output[ii][jj][kk]=input[kk*N2*N1 + jj*N1 + ii];
	// if(mympi.ThisTask == 0) printf("%d %d %d %f\n",ii,jj,kk,output[ii][jj][kk]);
      }
}

/**
 * Convert 3D-fftw_real array to 1D-float array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 3D-fftw_real array output
 * @param output 1D-float array input
 * @param Nion Number of the choices of total reionizing photons/neutral baryon
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void pack_4d_array_mpi_transfer(fftw_real ****input, float *output, int Nnion, int N1, int N2, int N3) { 
  int ii,jj,kk,jk;
  for(jk=0;jk<Nnion;jk++)
    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++) {
	  output[jk*N1*N2*N3 + kk*N1*N2 + jj*N1 + ii] = input[jk][ii][jj][kk];
	}
}
/** 
 * Convert 1D-float array to 3D-fftw_real array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 1D-float array input
 * @param output 3D-fftw_real array output
 * @param Nion Number of the choices of total reionizing photons/neutral baryon
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void unpack_4d_array_mpi_transfer(float *input, fftw_real ****output, int Nnion,int N1, int N2, int N3) {
  int ii,jj,kk,jk;
  for(jk=0;jk<Nnion;jk++)
    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++)
	  output[jk][ii][jj][kk]=input[jk*N1*N2*N3 + kk*N1*N2 + jj*N1 + ii];
}


/** 
 * Read in density of matter in cubep3m format (Fortran binary)
 *  
 *
 * @param filename Filename (char* input) 
 * @param buffer_3d Output buffer[N1*N2*N3]
 * @param robar_p Average density of baryons (output)
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 * @param vomegam Omega_matter (input)
 * @param vomegab Omega_baryon (input)
 */
void read_density(char *filename, float *buffer_3d, double *robar_p, int N1, int N2, int N3, float vomegam, float vomegab) {  
  int ii;
  int n1,n2,n3;
  FILE *inp;
  // printf("start read_density\n");
  if((inp=fopen(filename,"rb")) == NULL) {
    debug_checkpoint();
    printf("Cannot open density file: %s\nTerminating....\n",filename);
    exit(1);
  }
  *robar_p=0.;
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
  if(n1 != N1 || n2 != N2 || n3 != N3) {
    printf("Grid dimensions in %s are not the same as in config file\n",filename);
    printf("Density file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
    printf("Terminate\n");
    exit(1);
  }
  fread(buffer_3d,sizeof(float),n1*n2*n3,inp);
  fclose(inp);
  
  for(ii=0;ii<n1*n2*n3;ii++) {
    buffer_3d[ii] *= vomegab/vomegam;
    *robar_p += buffer_3d[ii];
  }
  *robar_p /= (1.*(n1)*(n2)*(n3));
}
void write_xfrac(char *dirname, char *z_out, float *buffer_4d, fftw_real ***nh, float robar, float *nion, int Nnion, int N1, int N2, int N3) {
  FILE *inp,*sinp;
  float xh1;
  int ii,jj,kk,jk,ll,start_ll;
  double t_start, t_stop;
  double vion[Nnion],roion[Nnion];
  char filename[2000];
  char buff[2000];
  int dummy;
  sinp = fopen(input_param.summary_file,"a");
  fprintf(sinp,"%s\t",z_out);
  for(jk=0;jk<Nnion;jk++) {
      t_start = Get_Current_time();
    
      vion[jk]=0.0;
      roion[jk]=0.0;
      // Defining the ionization map output file name
      // This is based on the value of nion assigned to it
      sprintf(buff,"mkdir -p %s/%4.2f",dirname,nion[jk]);
      system(buff);
      sprintf(filename,XFRACFILEPATTERN,dirname,nion[jk],PREFIX,z_out);
      printf("Saving %s\n",filename);
      ii=0; jj=0; kk=0;
      start_ll = jk*N1*N2*N3;
      for(ll=start_ll; ll<(jk+1)*N1*N2*N3; ll++) {

	xh1 = min(1.0,buffer_4d[ll]);
	// xh1 = max(xh1,0.0);
	buffer_4d[ll] = xh1;

	vion[jk] += xh1;
	roion[jk] += xh1*nh[ii][jj][kk];
	ii = (ll-start_ll) % N1;
	jj = ((ll-start_ll)/N1) % N2;
	kk = ((ll-start_ll)/(N1*N2)) % N3;
      }
      
      inp=fopen(filename,"wb+");
      // Writing the x_HI map in binary
      // In the begining 3 integers are written which defines the size
      // of the x_HI array
      dummy = sizeof(int)*3;
      fwrite(&dummy,sizeof(int),1,inp);
      fwrite(&N1,sizeof(int),1,inp);
      fwrite(&N2,sizeof(int),1,inp);
      fwrite(&N3,sizeof(int),1,inp);
      fwrite(&dummy,sizeof(int),1,inp);
      dummy = sizeof(float)*N1*N2*N3;
      fwrite(&dummy,sizeof(int),1,inp);
      fwrite(&buffer_4d[jk*N1*N2*N3],sizeof(float),N1*N2*N3,inp);
      fwrite(&dummy,sizeof(int),1,inp);
      fclose(inp);
      roion[jk]/=robar*(N1*N2*N3); // mass avg xHI
      vion[jk]/=(1.*N1*N2*N3); // volume avg xHI
      // roion[jk]/=(float)robar; // divide by H density to get mass avg. xHI
      t_stop = Get_Current_time();
      printf("nion = %f obtained vol. avg. x_HI=%e mass avg. x_HI=%e : %lf s\n",nion[jk],vion[jk],roion[jk],t_stop-t_start);
      fprintf(sinp,"%e\t%e\t",vion[jk],roion[jk]);
    }
  fprintf(sinp,"\n");
  fclose(sinp);
}

/** 
 * 
 * Read in the history of ionized fraction (previous snapshot).
 * 
 * @param dirname Directory of the xh (input) 
 * @param buffer_4d 4D float output
 * @param nion_list List of minimum N_gamma/N_b (input) 
 * @param Nnion Number of elements in nion_list
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 */
void read_xfrac(char *dirname, char *z, float *buffer_4d, float *nion, int Nnion, int N1, int N2, int N3) {
  FILE *inp;
  int n1,n2,n3;
  char filename[2048];
  int i,jk;
  float zval;
  int dummy;
  sscanf(z,"%f",&zval);
  if(zval < 0.0) {
    for(i=0;i<Nnion*N1*N2*N3;i++)
      buffer_4d[i] = 0.0;
    return;
  }
  for(jk=0;jk<Nnion;jk++) {
    sprintf(filename,XFRACFILEPATTERN,dirname,nion[jk],PREFIX,z);
    // printf("Reading %s\n",filename);
    if((inp=fopen(filename,"rb")) == NULL) {
      debug_checkpoint();
      printf("Cannot open %s\nTerminate...\n",filename);
      exit(1);
    }
    fread(&dummy,sizeof(int),1,inp);
    fread(&n1,sizeof(int),1,inp);
    fread(&n2,sizeof(int),1,inp);
    fread(&n3,sizeof(int),1,inp);
    fread(&dummy,sizeof(int),1,inp);
    if(n1 != N1 || n2 != N2 || n3 != N3) {
      printf("Grid dimensions in %s are not the same as in config file\n",filename);
      printf("Xfrac file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
      printf("Terminate\n");
      exit(1);
    }
    fread(&dummy,sizeof(int),1,inp);
    fread(&buffer_4d[jk*N1*N2*N3],sizeof(float),n1*n2*n3,inp);
    fread(&dummy,sizeof(int),1,inp);
    fclose(inp);
  }
}


/** 
 * Read sources list. The data will be the SFR in cubep3m_gridmass/year
 * 
 * @param filename Filename (char* input) 
 * @param buffer_3d Output buffer[N1*N2*N3]
 * @param robar_p Average density of photons (output)
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 */
void read_sources(char *filename, float *buffer_3d, double *robarhalo_p, int N1, int N2, int N3) {
  FILE *inp;
  int ii;
  int n1,n2,n3;
  
  if((inp=fopen(filename,"rb")) == NULL) {
    debug_checkpoint();
    printf("Cannot open sources file: %s\nTerminating....\n",filename);
    exit(1);
  }
  *robarhalo_p=0.;
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
  if(n1 != N1 || n2 != N2 || n3 != N3) {
    printf("Grid dimensions in %s are not the same as in config file\n",filename);
    printf("Source file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
    printf("Terminate\n");
    exit(1);
  }
  fread(buffer_3d,sizeof(float),n1*n2*n3,inp);
  fclose(inp);
  for(ii=0;ii<n1*n2*n3;ii++) {
    *robarhalo_p += buffer_3d[ii];
  }
  *robarhalo_p /= (1.*(n1)*(n2)*(n3));
}

