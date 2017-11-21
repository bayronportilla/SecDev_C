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

/*
double de_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);
*/

int main (void){

  FILE *fp;
  fp = fopen("data.dat","w");

  Inpar st;
  st = params();

  double mu      = 1.0;

  
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 16);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (16);
  gsl_odeiv2_system sys = { func, NULL, 16, &mu}; // Define sistema de ecuaciones
  
  double h = 1e-6;
  
  Inpar stc;
  stc = ConverToCan(st);

  
  double y[16] = { stc.a_in, stc.a_out, stc.e_in, stc.e_out, 
		   stc.I_in, stc.I_out, stc.W_in, stc.W_out,
		   stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
		   stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz}; // y[number of entries of the array] = {}
  
  /*
  printf("%e \n",dW_in_dt(stc.a_in,stc.a_out,stc.e_in,
			  stc.e_out,stc.I_in,stc.I_out,
			  stc.W_in,stc.W_out,stc.w_in,
			  stc.w_out,stc.Om_Ax,stc.Om_Ay,
			  stc.Om_Az,stc.Om_Bx,stc.Om_By,
			  stc.Om_Bz,0,stc));
    
  printf("%e \n",dW_in_dt(st.a_in,st.a_out,st.e_in,
			  st.e_out,st.I_in,st.I_out,
			  st.W_in,st.W_out,st.w_in,
			  st.w_out,st.Om_Ax,st.Om_Ay,
			  st.Om_Az,st.Om_Bx,st.Om_By,
  			  st.Om_Bz,0,st));

  */

  double t = stc.t_ini;
  while (t < stc.t_end)
    {
      int status = gsl_odeiv2_evolve_apply (e,c,s,&sys,&t,stc.t_end,&h,y);

      if (status != GSL_SUCCESS)
          break;

      printf ("%.5e %.5e %.5e\n", t, y[0], y[15]);
    }
  
  

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  fclose(fp);


  
  return 0;
  
}

