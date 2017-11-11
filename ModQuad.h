int func (double t, const double y[], double f[], void *params){
  (void)(t); // avoid unused parameter warning 
  double mu = *(double *)params;
  double x  = y[0];
  double vx = y[1];
  f[0] = vx;
  f[1] = -mu*mu*x;
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
