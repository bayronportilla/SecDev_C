#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

int main(void){

  ////////////////////////////////////////////////////////////
  //
  // vectors
  //
  ////////////////////////////////////////////////////////////
  
  gsl_vector *v = gsl_vector_alloc (3);

  gsl_vector_set (v, 0, 0.23);
  gsl_vector_set (v, 1, 1.23);
  gsl_vector_set (v, 2, 2.23);
  
  printf ("v_%d = %g\n", 0, gsl_vector_get (v, 0));
  printf ("v_%d = %g\n", 1, gsl_vector_get (v, 1));
  printf ("v_%d = %g\n", 2, gsl_vector_get (v, 2));
  
  ////////////////////////////////////////////////////////////
  //
  // matrices
  //
  ////////////////////////////////////////////////////////////
  
  int i, j;
  gsl_matrix *m = gsl_matrix_alloc (3, 3);

  gsl_matrix_set (m, 0, 0, 3);
  gsl_matrix_set (m, 0, 1, 0);
  gsl_matrix_set (m, 0, 2, 1);
  gsl_matrix_set (m, 1, 0, 2);
  gsl_matrix_set (m, 1, 1, 1);
  gsl_matrix_set (m, 1, 2, 0);
  gsl_matrix_set (m, 2, 0, 10);
  gsl_matrix_set (m, 2, 1, 9);
  gsl_matrix_set (m, 2, 2, 7);

  
  ////////////////////////////////////////////////////////////
  //
  // multiplication
  //
  ////////////////////////////////////////////////////////////
  
  gsl_vector *y = gsl_vector_alloc (3); // Result will be stored in y

  gsl_blas_dgemv(CblasTrans, 1.0, m, v, 0.0, y);

  
  printf ("v_%d = %f\n", 0, gsl_vector_get (y, 0));
  printf ("v_%d = %f\n", 1, gsl_vector_get (y, 1));
  printf ("v_%d = %f\n", 2, gsl_vector_get (y, 2));

  
  gsl_vector_free (v);
  gsl_vector_free (y);
  gsl_matrix_free (m);
  
  
  return 0;


}
