#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main (void){
  
  FILE *file;
  file = fopen("radius.dat","r");

  int j = 0;
  int Nlines = 0;
  float col1,col2;
  double *tim,*rad;
  
  tim = (double *)malloc(1*sizeof(double));
  rad = (double *)malloc(1*sizeof(double));

  while(fscanf(file,"%f %f",&col1,&col2)!=EOF){
    tim[j] = col1;
    rad[j] = col2;
    j++;
    tim = (double *)realloc(tim,sizeof(double)*(j+1));
    rad = (double *)realloc(rad,sizeof(double)*(j+1));
  }

  Nlines = j;

  
  ////////////////////////////////////////////////////////////
  // Start interpolation
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, Nlines);
  gsl_spline_init(spline, tim, rad, Nlines);
    
  double ti,yi;

  for (ti = tim[0]; ti < tim[Nlines-1]; ti += 0.01){
    yi = gsl_spline_eval (spline, ti, acc);
    printf ("%g %g\n", ti, yi);
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(tim);
  free(rad);  

  return 0;
}
