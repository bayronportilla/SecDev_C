#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define PI 3.14159
  
int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -mu*mu*y[0];
  return GSL_SUCCESS;
}



int jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
  (void)(t); // avoid unused parameter warning 
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}



int main (void){

  FILE *fp;
  fp = fopen("data.dat","w");
  double mu = 1.0;
  gsl_odeiv2_system sys = { func, jac,2, &mu };

  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys,gsl_odeiv2_step_rk4,1e-3,1e-8,1e-8);  //Driver system

  double t_ini   = 0.0; //t_ini
  double t_end   = 5*2*PI/mu;
  double n_steps = 1000;
  double h       = (t_end-t_ini)/n_steps;
  double y[2]    = { 1.0, 0.0 }; // y[number of entries of the array] = {}
  int i, s;

  double t = t_ini;

  /*
  printf("time = %f\n",t_ini);
  s = gsl_odeiv2_driver_apply_fixed_step (d, &t_ini, h, 1, y); // driver system, t_ini, step size, n_steps, state vector
  printf("time = %f\n",t_ini);
  */
    
  
  for (i = 0; i < n_steps; i++){
    
    s = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y); // driver system, t_ini, step size, n_steps, state vector

    if (s != GSL_SUCCESS){
      printf ("error: driver returned %d\n", s);
      break;
    }

    //printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    fprintf(fp,"%.5e %.5e %.5e\n",t, y[0], y[1]);
  }
  
  fclose(fp);

  gsl_odeiv2_driver_free (d);
  return s;

  
}

