/**
 * @file   ionz_funcs.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sun Oct 12 00:59:24 2014
 * 
 * @brief  Several core functions
 * 
 * 
 */

#include "ion.h"

/** 
 * Set up memory for nh, ngamma and nxion
 * 
 * @param Nnion Number of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void Setting_Up_Memory_For_ionz(int Nnion, int N1, int N2, int N3) {

}

/** 
 * Smooth radius with FFT
 * 
 * @param ro_dum density field
 * @param Radii Smooth radius 
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void smooth(fftw_real ***ro_dum,float Radii,int N1,int N2, int N3) {
  int i,j,k,index;
  float tempre,tempim,tot;
  fftw_real ***rosp;
  fftw_complex *A;
  fftw_complex *B;

  rfftwnd_plan p_ro; // for FFT
  rfftwnd_plan q_ro; // for FFT
  
  rosp = allocate_fftw_real_3d(N1,N2,N3+2);
  /****************************************************/	
  /* Creating the plans for forward and reverse FFT's */
  
  p_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);  
  q_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);
 
  //generating the filtering function
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
  	rosp[i][j][k]=0.0;
 
  //Radii is radius of the sphere in grid unit
  //generating a sphere at the centre of the box
 
  tot=0.;
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++) {
	if((float)((N1/2-i)*(N1/2-i)+(N2/2-j)*(N2/2-j)+(N3/2-k)*(N3/2-k))<=Radii*Radii)
	  rosp[i][j][k]=1.0;//centre N1/2,N2/2,N3/2
	tot += rosp[i][j][k];
      }
  //Sphere generation complete 
  //Doing Fourier Transform of the sphere
  rfftwnd_one_real_to_complex(p_ro,&rosp[0][0][0], NULL);
  B=(fftw_complex*)&(rosp[0][0][0]);

  //We will multiply the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to one corner of the box from box centre while applying boundary condition below
  //----------------------------------------------------------------------

  //Doing Fourier Transform of the density field
  rfftwnd_one_real_to_complex(p_ro,&ro_dum[0][0][0], NULL);
  A=(fftw_complex*)&(ro_dum[0][0][0]);
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<=N3/2;k++)    { 
	index = i*N2*(N3/2 +1) + j*(N3/2 +1) + k;
	tempre=(A[index].re*B[index].re-A[index].im*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	tempim=(A[index].im*B[index].re+A[index].re*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	//multiplying the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to box centre from one corner of the box after complex to real FT
	A[index].re=tempre;
	A[index].im=tempim;
      }
  rfftwnd_one_complex_to_real(q_ro,(fftw_complex *) &ro_dum[0][0][0], NULL);
  free_fftw_real_3d(rosp,N1,N2,N3+2);
  rfftwnd_destroy_plan(p_ro);
  rfftwnd_destroy_plan(q_ro);
  /* A and B are aliases so there is no need to free them... Boyd */
  // fftw_free(A);
  // fftw_free(B);
}


void subgrid_reionization(fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, double robar, float *nion_p, int Nnion, int N1, int N2, int N3) {
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  int ii,jj,kk;
#else
  int len=N1*N2*(N3+2);
#endif
  int jk;
  for(jk=0;jk<Nnion;jk++) {
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
    //calculating avg. ionization frction
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	    nxion_p[jk][ii][jj][kk]=min(nion_p[jk]*ngamma_p[ii][jj][kk]/nh_p[ii][jj][kk],1.0);
	}
#else
    fortran_subgrid_reionization(&nh_p[0][0][0],&ngamma_p[0][0][0],&nxion_p[jk][0][0][0],&nion_p[jk],&len);    
#endif
  }
}

/** 
 * Do subgrid semi-numnerical reionization estimation
 * 
 * @param Radii Radius 
 * @param nh_p Pointer to baryon density field
 * @param ngamma_p Pointer to photon density field
 * @param xfrac_p Pointer to old Xfrac field
 * @param nxion_p Pointer to Xfrac field
 * @param nion_p nion arrays
 * @param Nnion_p Total element of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void subgrid_reionization_with_xfrac(fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****xfrac_p, fftw_real ****nxion_p, double robar, float *nion_p, int Nnion, int N1, int N2, int N3) {
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  int ii,jj,kk;
  double nh;
#else
  int len = N1*N2*(N3+2);
#endif
  int jk;
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  for(jk=0;jk<Nnion;jk++) {    
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	  nh = nh_p[ii][jj][kk]*(1.0-xfrac_p[jk][ii][jj][kk]);
	  nxion_p[jk][ii][jj][kk]=min(1.0,xfrac_p[jk][ii][jj][kk] + (1.0-xfrac_p[jk][ii][jj][kk])*nion_p[jk]*ngamma_p[ii][jj][kk]/nh);	  
	}
  }
#else
  for(jk=0;jk<Nnion;jk++) {
    fortran_subgrid_reionization_with_xfrac(&nh_p[0][0][0],&ngamma_p[0][0][0],&xfrac_p[jk][0][0][0],&nxion_p[jk][0][0][0],&nion_p[jk],&len);
  }
#endif
}


/** 
 * Core function to do semi-numnerical reionization estimation
 * 
 * @param Radii Radius 
 * @param nh_p Pointer to baryon density field
 * @param ngamma_p Pointer to photon density field
 * @param xfrac_p Pointer to old Xfrac field
 * @param nxion_p Pointer to Xfrac field
 * @param nion_p nion arrays
 * @param Nnion Total element of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void reionization_with_xfrac(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****xfrac_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3) {
  fftw_real ***nhs,***ngammas;
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  int ii,jj,kk;
#else
  int len=N1*N2*(N3+2);
  float one = 1.;
#endif
  int jk;
  
  nhs=allocate_fftw_real_3d(N1,N2,N3+2);
  ngammas=allocate_fftw_real_3d(N1,N2,N3+2);

  for(jk=0;jk<Nnion;++jk) { 
    //Filling smoothing arrays with the dark matter and source density data
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	  nhs[ii][jj][kk]=nh_p[ii][jj][kk]*(1.0-xfrac_p[jk][ii][jj][kk]);
	  ngammas[ii][jj][kk]=ngamma_p[ii][jj][kk];	     
	}
#else
    fortran_prepar_fftw_real_3d_with_xfrac(&nh_p[0][0][0],&nhs[0][0][0],&xfrac_p[jk][0][0][0],&len);
    fortran_multiply_constant_fftw_real_3d(&ngamma_p[0][0][0],&one,&ngammas[0][0][0], &len);
#endif
    //Smoothing with real space spherical filter
    smooth(nhs,Radii,N1,N2,N3);
    smooth(ngammas,Radii,N1,N2,N3); 

#ifndef USE_FORTRAN_SPEEDUP_ARRAY
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	  //Checking the ionization condition
	  if(nhs[ii][jj][kk]<ngammas[ii][jj][kk]*nion_p[jk]) {
	    nxion_p[jk][ii][jj][kk]=1.;
	  }
	}
#else
    fortran_condition_ionize(&nhs[0][0][0],&ngammas[0][0][0],&nion_p[jk],&nxion_p[jk][0][0][0],&len);
#endif
  }
  free_fftw_real_3d(nhs,N1,N2,N3+2);
  free_fftw_real_3d(ngammas,N1,N2,N3+2);
}

/** 
 * Core function to do semi-numnerical reionization estimation (with Xfrac)
 * 
 * @param Radii Radius 
 * @param nh_p Pointer to baryon density field
 * @param ngamma_p Pointer to photon density field
 * @param nxion_p Pointer to Xfrac field
 * @param nion_p nion arrays
 * @param Nnion Total element of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3) {
  fftw_real ***nhs,***ngammas;
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  int ii,jj,kk;
#else
  int len = N1*N2*(N3+2);
  float one = 1.;
#endif
  int jk;

  nhs=allocate_fftw_real_3d(N1,N2,N3+2);
  ngammas=allocate_fftw_real_3d(N1,N2,N3+2);
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  for(ii=0;ii<N1;ii++)
    for(jj=0;jj<N2;jj++)
      for(kk=0;kk<N3;kk++) {
	//Filling smoothing arrays with the dark matter and source density data
	nhs[ii][jj][kk]=nh_p[ii][jj][kk];
	ngammas[ii][jj][kk]=ngamma_p[ii][jj][kk];
      }
#else
    fortran_multiply_constant_fftw_real_3d(&ngamma_p[0][0][0],&one,&ngammas[0][0][0], &len);
    fortran_multiply_constant_fftw_real_3d(&nh_p[0][0][0],&one,&nhs[0][0][0], &len);
#endif
  // printf("starting smoothing for radius of size %e (in units of grid size)\n",Radii);
  
  //Smoothing with real space spherical filter
  
  smooth(nhs,Radii,N1,N2,N3);
  smooth(ngammas,Radii,N1,N2,N3);
 
#ifndef USE_FORTRAN_SPEEDUP_ARRAY
  for(jk=0;jk<Nnion;jk++)
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	  if(nhs[ii][jj][kk]<nion_p[jk]*ngammas[ii][jj][kk])
	    nxion_p[jk][ii][jj][kk]=1.;
	}
#else
  for(jk=0;jk<Nnion;jk++)
    fortran_condition_ionize(&nhs[0][0][0],&ngammas[0][0][0],&nion_p[jk],&nxion_p[jk][0][0][0],&len);
#endif

  free_fftw_real_3d(nhs,N1,N2,N3+2);
  free_fftw_real_3d(ngammas,N1,N2,N3+2);
}




