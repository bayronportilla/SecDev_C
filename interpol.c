#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
int main (void){

  FILE *file;
  file = fopen("radius.dat","r");
  
  int j = 0;
  float a,b;
  int size_array = 99;
  double rad[size_array];
  double tim[size_array];

  if (file==NULL){
    printf("Error: can't open file.\n");
    return 1;
  }
  
  else{
    while(fscanf(file,"%f %f",&a,&b)!=EOF){
      tim[j] = a;
      rad[j] = b;
      j += 1;
    }
  }

  int k;
  /*
  for(k=0;k<size_array;k++){
    printf("%f\n",tim[k]);
  }
  */
  /*
  int l;
  FILE *fp2;
  fp2 = fopen("radio_interpol.dat","w");
  for(l=0;l<j;l++){
    fprintf(fp2,"%f   %f \n",tim[l],rad[l]);
  }
  */
  
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, size_array);
  gsl_spline_init(spline, tim, rad, size_array);
  
  double ti,yi;


  for (ti = tim[0]; ti < tim[size_array-1]; ti += 0.1){
    yi = gsl_spline_eval (spline, ti, acc);
    printf ("%g %g\n", ti, yi);

  }


  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  
  
  /*
  int i;
  double xi, yi, x[10], y[10];

  for (i = 0; i < 10; i++)
    {
      x[i] = i + 0.5 * sin (i);
      y[i] = i + cos (i * i);
    }

  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    //gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, 10);
    gsl_spline_init (spline, x, y, 10);

    for (xi = x[0]; xi < x[9]; xi += 0.5)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        printf ("%g %g\n", xi, yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
  */
  return 0;
}
