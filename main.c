#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Units.h"
#include "Params.h"
#include "ModQuad.h"
#include "ConverCan.h"
#include <libconfig.h>

#define PI 3.14159


int main (void){

  FILE *fp;
  fp = fopen("data.dat","w");

  Inpar st;
  st = params();
  Inpar stc;
  stc = ConverToCan(st);

  
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
  gsl_odeiv2_system sys = { func, NULL, 16, &st}; // Define sistema de ecuaciones
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-3, 1e-8, 1e-8);
  
  const double h = 1.0e3;


  
  double y[16] = { stc.a_in, stc.a_out, stc.e_in, stc.e_out, 
		   stc.I_in, stc.I_out, stc.W_in, stc.W_out,
		   stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
		   stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz}; // y[number of entries of the array] = {}
  

  
  int i, s;
  double t = stc.t_ini;

  /*
  for (i = 0; i < 100; i++){
    s = gsl_odeiv2_driver_apply_fixed_step (d, &t, 1e3, 10, y);
    
    if (s != GSL_SUCCESS)
      {
	printf ("error: driver returned %d\n", s);
	break;
      }
    
    printf("%.5e %.5e %.5e\n", t, y[2], y[4]+y[5]);
    fprintf(fp,"%.5e %.5e %.5e\n", t, y[2], y[4]+y[5]);
  }
  */



  while(t<stc.t_end){
    s = gsl_odeiv2_driver_apply_fixed_step (d, &t, 1e3, 10, y);
    
    if (s != GSL_SUCCESS)
      {
	printf ("error: driver returned %d\n", s);
	break;
      }
    
    printf("%.5e %.5e %.5e\n", t, y[2], y[4]+y[5]);
    fprintf(fp,"%.5e %.5e %.5e\n", t, y[2], y[4]+y[5]);
  }

  



  gsl_odeiv2_driver_free (d);

  fclose(fp);


  
  return 0;
  
}

