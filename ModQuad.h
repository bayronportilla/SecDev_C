////////////////////////////////////////////////////////////
//
// This matrix takes me from the inertial system to the
// orbital system
//
////////////////////////////////////////////////////////////

double *RotMat(double W, double I, double w, double *in_array, char dir[]){

  /*
    Esta funcion permite hacer la conversion de un arreglo unidimensional de 
    tres entradas medidas respecto a un sistema de referencia inicial a otro 
    arreglo de la misma dimension con entradas medidas en otro sistema que se
    ha construido mediante tres rotaciones sucesivas parametrizadas por los 
    angulos W,I,w. Tambien permite hacer la conversion inversa. En el contexto 
    de este programa, las posibles direcciones son: de un sistema inercial a uno 
    orbital "InOrb" y de uno orbital a uno inercial "OrbIn". Para usarlo, los
    argumentos a introducir (respetando el orden son) W, I, w, un puntero a un vector
    que contiene las entradas en el plano de partida (inercial u orbital) y una cadena
    de caracteres que indica la direccion de la conversion "InOrb" o "OrbIn". El return 
    de esta funcion es un vector unidimensional de tres entradas, medidas en el plano de 
    llegada. 

    Ej 1. 

    #include <ModQuad.h>
    ... y demas includes de GSL y math.h

    int main(void){
       double vector[3]={0.1,0.2,0.3};
       printf("%f n",RotMat(2,3,4,vector,"InOrb")[0]);
    }
    
   */
  
  ////////////////////////////////////////////////////////////
  //
  // vectors
  //
  ////////////////////////////////////////////////////////////

  static double out_array[3]; // the output array with the elements in the orbital plane
  
  gsl_vector *v = gsl_vector_alloc (3);

  gsl_vector_set (v, 0, in_array[0]);
  gsl_vector_set (v, 1, in_array[1]);
  gsl_vector_set (v, 2, in_array[2]);
  
  ////////////////////////////////////////////////////////////
  //
  // matrices
  //
  ////////////////////////////////////////////////////////////
  
  gsl_matrix *m = gsl_matrix_alloc (3, 3);
  
  gsl_matrix_set (m, 0, 0,  cos(W)*cos(w) - sin(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 1,  sin(W)*cos(w) + cos(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 2,  sin(I)*sin(w));
  gsl_matrix_set (m, 1, 0, -cos(W)*sin(w) - sin(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 1, -sin(W)*sin(w) + cos(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 2,  sin(I)*cos(w));
  gsl_matrix_set (m, 2, 0,  sin(I)*sin(W));
  gsl_matrix_set (m, 2, 1, -sin(I)*cos(W));
  gsl_matrix_set (m, 2, 2,  cos(I));
 
  ////////////////////////////////////////////////////////////
  //
  // multiplication
  //
  ////////////////////////////////////////////////////////////
  
  gsl_vector *y = gsl_vector_alloc (3); // Result will be stored in y

  char str1[] = "InOrb";
  char str2[] = "OrbIn";
  
  if( strcmp(dir,str1) == 0){
      gsl_blas_dgemv(CblasNoTrans, 1.0, m, v, 0.0, y);
    }

  if( strcmp(dir,str2) == 0){
      gsl_blas_dgemv(CblasTrans, 1.0, m, v, 0.0, y);
    }
  
  out_array[0] = gsl_vector_get (y, 0);
  out_array[1] = gsl_vector_get (y, 1);
  out_array[2] = gsl_vector_get (y, 2);

  gsl_vector_free (v);
  gsl_vector_free (y);
  gsl_matrix_free (m);
 
  return out_array;
  
}



////////////////////////////////////////////////////////////
//
// Definition of generic functions
//
////////////////////////////////////////////////////////////

double f_1(double e){
  return ( 1.0 + 3.0*pow(e,2) + (3.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}

double f_2(double e){
  return ( 1.0 + (15.0/2.0)*pow(e,2) + (45.0/8.0)*pow(e,4) + (5.0/16.0)*pow(e,6) ) / pow(1.0-e*e,13.0/2.0);
}

double f_3(double e){
  return ( 1.0 + (9.0/2.0)*pow(e,2) + (5.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}

double f_4(double e){
  return ( 1.0 + (3.0/2.0)*pow(e,2) + (1.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}

double f_5(double e){
  return ( 1.0 + (15.0/4.0)*pow(e,2) + (15.0/8.0)*pow(e,4) + (5.0/64.0)*pow(e,6) ) / pow(1.0-e*e,13.0/2.0);
}

double G = 6.67e-11;

double n(double m_1, double m_2, double a){
  return pow(G*(m_1+m_2)/pow(a,3),0.5);
}

double L_1(double m_1, double m_2, double a_1){
  return m_1*m_2/(m_1+m_2) * pow(G*(m_1+m_2)*a_1,0.5);
}


double L_2(double m_1,double m_2,double m_3,double a_2){
  return m_3*(m_1+m_2)/(m_1+m_2+m_3) * pow(G*(m_1+m_2+m_3)*a_2,0.5);
}


double G_1(double m_1,double m_2,double a_1,double e_1){
  return L_1(m_1,m_2,a_1)*pow(1.0-e_1*e_1,0.5);
}

double G_2(double m_1,double m_2,double m_3,double a_2,double e_2){
  return L_2(m_1,m_2,m_3,a_2)*pow(1.0-e_2*e_2,0.5);
}

double H_1(double m_1,double m_2,double a_1,double e_1,double I_1){
  return G_1(m_1,m_2,a_1,e_1)*cos(I_1);
}

double H_2(double m_1,double m_2,double m_3,double a_2,double e_2,double I_2 ){
  return G_2(m_1,m_2,m_3,a_2,e_2)*cos(I_2);
}

double C2(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2){
  return G*G/16.0 * pow(m_1+m_2,7)/pow(m_1+m_2+m_3,3) * pow(m_3,7)/pow(m_1*m_2,3) * pow(L_1(m_1,m_2,a_1),4)/pow(L_2(m_1,m_2,m_3,a_2),3) * pow(G_2(m_1,m_2,m_3,a_2,e_2),3);
}

double epsilon_m(double m_1,double m_2,double a_1,double a_2,double e_2){
  return (m_1-m_2)/(m_1+m_2) * a_1/a_2 * e_2/(1.0-e_2*e_2);
}

double C3(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2){
  return -C2(m_1,m_2,m_3,a_1,a_2,e_2)*15.0/4.0 * epsilon_m(m_1,m_2,a_1,a_2,e_2)/e_2;
}

double tf_A(double tv_A, double R_A, double k_A,double a,double m_A, double m_B){
  return tv_A/9.0 * pow(a/R_A,8) * m_A*m_A/((m_A+m_B)*m_B) * 1.0/pow(1.0+2.0*k_A,2);
}

double tf_B(double tv_B, double R_B, double k_B, double a, double m_A, double m_B){
  return tv_B/9.0 * pow(a/R_B,8) * m_B*m_B/((m_A+m_B)*m_A) * 1.0/pow(1.0+2.0*k_B,2);
}

double V_A(double tv_A, double R_A, double k_A, double m_A, double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az){
  
  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  
  return 9.0/tf_A(tv_A,R_A,k_A,a,m_A,m_B) * (f_5(e) - 11.0/18.0*Om_Az_orb/n(m_A,m_B,a) * f_4(e));
}

double V_B(double tv_B, double R_B, double k_B, double m_A, double m_B, 
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz){
  
  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  
  return 9.0/tf_B(tv_B,R_B,k_B,a,m_A,m_B) * (f_5(e) - 11.0/18.0*Om_Bz_orb/n(m_A,m_B,a) * f_4(e));
}


double W_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){
  
  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];

  return 1.0/tf_A(tv_A,R_A,k_A,a,m_A,m_B) * ( f_2(e) - Om_Az_orb/n(m_A,m_B,a) * f_1(e) );
}


double W_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){
  
  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];

  return 1.0/tf_B(tv_B,R_B,k_B,a,m_A,m_B) * ( f_2(e) - Om_Bz_orb/n(m_A,m_B,a) * f_1(e) );
}


double X_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){

  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
  return -m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Az_orb*Om_Ax_orb/pow(1.0-e*e,2) - Om_Ay_orb/(2.0*n(m_A,m_B,a)*tf_A(tv_A,R_A,k_A,a,m_A,m_B)) * f_3(e);
}

double X_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
  return -m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Bz_orb*Om_Bx_orb/pow(1.0-e*e,2) - Om_By_orb/(2.0*n(m_A,m_B,a)*tf_B(tv_B,R_B,k_B,a,m_A,m_B)) * f_3(e);
}


double Y_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){

  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
 
  return -m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Az_orb*Om_Ay_orb/pow(1.0-e*e,2) + Om_Ax_orb/(2.0*n(m_A,m_B,a)*tf_A(tv_A,R_A,k_A,a,m_A,m_B)) * f_4(e);
  
}

double Y_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
 
  return -m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Bz_orb*Om_By_orb/pow(1.0-e*e,2) + Om_Bx_orb/(2.0*n(m_A,m_B,a)*tf_B(tv_B,R_B,k_B,a,m_A,m_B)) * f_4(e);
  
}

double Z_A(double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){

  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
  return m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * ( (2.0*pow(Om_Az_orb,2) - pow(Om_Ay_orb,2) - pow(Om_Ax_orb,2))/(2.0*pow(1.0-e*e,2)) + 15.0*G*m_B/pow(a,3) * f_4(e) );
  
}


double Z_B(double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
  return m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * ( (2.0*pow(Om_Bz_orb,2) - pow(Om_By_orb,2) - pow(Om_Bx_orb,2))/(2.0*pow(1.0-e*e,2)) + 15.0*G*m_A/pow(a,3) * f_4(e) );
  
}




////////////////////////////////////////////////////////////
//
// Differential equations
//
////////////////////////////////////////////////////////////


double da_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){

  double da_in_dt_orb;
  double da_in_dt_tid;

  da_in_dt_orb = 0.0;
  
  /*
  da_in_dt_tid = -2.0*a_in*( W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)
	         + W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) - \
                 2.0*a_in*e_in**2/(1.0-e_in**2) * ( V_A(tv_A,R_A,m_A,m_B,k_A,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + \
	         V_B(tv_B,R_B,m_A,m_B,k_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) );
  */
  
  //da_in_dt_tid = W_A(a_in);
  return 0;  
}


/*
double da_in_dt(Inpar params){
  
  double da_in_dt_orb = 0.0;
  double da_in_dt_tid = 0.0;
  
  return (params.q_orb*2.0) + (params.q_tid*3.0);
}
*/

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
