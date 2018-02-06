#include "allvars.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <string.h>

double interpol(Inpar st, double t, char quantity[]){

  FILE *file;
  file = fopen("radius.dat","r");

  int j = 0;
  int Nlines = 0;
  float col1,col2,col3,col4;
  double *tim_A,*tim_B,*rad_A,*rad_B;
  
  tim_A = (double *)malloc(1*sizeof(double));
  tim_B = (double *)malloc(1*sizeof(double));
  rad_A = (double *)malloc(1*sizeof(double));
  rad_B = (double *)malloc(1*sizeof(double));

  while(fscanf(file,"%f %f %f %f",&col1,&col2,&col3,&col4)!=EOF){
    tim_A[j] = col1;
    tim_B[j] = col2;
    rad_A[j] = col3;
    rad_B[j] = col4;
    j++;
    tim_A = (double *)realloc(tim_A,sizeof(double)*(j+1));
    tim_B = (double *)realloc(tim_B,sizeof(double)*(j+1));
    rad_A = (double *)realloc(rad_A,sizeof(double)*(j+1));
    rad_B = (double *)realloc(rad_B,sizeof(double)*(j+1));
  }

  Nlines = j;
  
  ////////////////////////////////////////////////////////////
  //
  // Finding max and min of time entries
  //
  ////////////////////////////////////////////////////////////

  int i;
  gsl_vector *time_A_gsl = gsl_vector_alloc (Nlines);
  gsl_vector *time_B_gsl = gsl_vector_alloc (Nlines);
  gsl_vector *time_compa = gsl_vector_alloc (2);
  
  for (i = 0; i < Nlines; i++){
      gsl_vector_set (time_A_gsl, i, tim_A[i]);
      gsl_vector_set (time_B_gsl, i, tim_B[i]);
  }

  double min_time_A, min_time_B;
  double max_time_A, max_time_B;

  min_time_A = gsl_vector_min(time_A_gsl);
  min_time_B = gsl_vector_min(time_B_gsl);
  max_time_A = gsl_vector_max(time_A_gsl);
  max_time_B = gsl_vector_max(time_B_gsl);
  
  if (min_time_A<min_time_B){
    gsl_vector_set (time_compa,0,min_time_A);
  }
  else{
    gsl_vector_set (time_compa,0,min_time_B);
  }

  if (max_time_A>max_time_B){
    gsl_vector_set (time_compa,1,max_time_A);
  }
  else{
    gsl_vector_set (time_compa,1,max_time_B);
  }


  
  ////////////////////////////////////////////////////////////
  //
  // Printing info
  //
  ////////////////////////////////////////////////////////////

  //printf("t_ini can't be lower than:   %1.3e yr\n",gsl_vector_get(time_compa,0)*1e9);
  //printf("t_end can't be greater than: %1.3e yr\n",gsl_vector_get(time_compa,1)*1e9);
  

  ////////////////////////////////////////////////////////////
  //
  // Start interpolation
  //
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Body A
  
  gsl_interp_accel *acc_rad_A = gsl_interp_accel_alloc ();
  gsl_spline *spline_rad_A    = gsl_spline_alloc (gsl_interp_linear, Nlines);
  gsl_spline_init(spline_rad_A,tim_A,rad_A,Nlines);
  
  if (strcmp(quantity,"radius_A")==0){
    return gsl_spline_eval (spline_rad_A, t, acc_rad_A);
  }


  
  ////////////////////////////////////////////////////////////
  // Body B
  
  gsl_interp_accel *acc_rad_B = gsl_interp_accel_alloc ();
  gsl_spline *spline_rad_B    = gsl_spline_alloc (gsl_interp_linear, Nlines);
  gsl_spline_init(spline_rad_B,tim_B,rad_B,Nlines);
  
  if (strcmp(quantity,"radius_B")==0){
    return gsl_spline_eval (spline_rad_B, t, acc_rad_B);
  }



  ////////////////////////////////////////////////////////////
  //
  // Releasing memory
  //
  ////////////////////////////////////////////////////////////
  
  gsl_spline_free (spline_rad_A);
  gsl_spline_free (spline_rad_B);
  gsl_interp_accel_free (acc_rad_A);
  gsl_interp_accel_free (acc_rad_B);
  free(tim_A);
  free(tim_B);
  free(rad_A);
  free(rad_B);  
  
  
  return 0;
  
}
