#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Units.h"
#include "Params.h"
#include "ModQuad.h"

#include <libconfig.h>

#define PI 3.14159


int main (void){

  FILE *fp;
  fp = fopen("data.dat","w");
  
  Inpar st;
  st = params();
  
  double mu      = 1.0;
  double t_ini   = st.t_ini;
  double t_end   = st.t_end;
  int    n_steps = st.n_steps;

  gsl_odeiv2_system sys = { func, jac, 2, &mu };
  gsl_odeiv2_driver *d  = gsl_odeiv2_driver_alloc_y_new (&sys,gsl_odeiv2_step_rk4,1e-3,1e-8,1e-8);  //Driver system

  double h = (t_end-t_ini)/n_steps;
  double y[2] = { 1.0, 0.0 }; // y[number of entries of the array] = {}

  int i, s;
  double t = t_ini;


  /*
  printf("%e \n",Z_B(1e6,0.25,1.989e30,5.97e24,
		     1e11,0.2,1.1,1.1,1.1,1e-5,1e-5,1e-5));
  */

  
  printf("%e \n",dOm_Bz_dt(st.a_in,st.a_out,st.e_in,
			   st.e_out,1.1292280260403313,0.005235987755982988,
			   st.W_in,st.W_out,st.w_in,
			   st.w_out,1e-5,1e-5,
			   1e-5,1e-5,1e-5,
			   1e-5,0.0,st));
  

  /*
  printf("%e \n",W_B(st.tv_A,st.R_A,st.k_A,
		     st.m_A,st.m_B,st.a_in,
		     st.e_in,st.W_in,2.2,
		     st.w_in,1e-5,1e-5,
		     1e-5));
  */
  
  /*
 
  for (i = 0; i < n_steps; i++){    
    s = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y); // driver system, t_ini, step size, n_steps, state vector
    if (s != GSL_SUCCESS){
      printf ("error: driver returned %d\n", s);
      break;
    }
    fprintf(fp,"%.5e %.5e %.5e\n",t, y[0], y[1]);
  }
  */
  

  fclose(fp);
  gsl_odeiv2_driver_free (d);

  return s;
  
}

