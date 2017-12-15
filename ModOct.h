

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


double G = 1.0;


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
  return G*G/16.0 * pow(m_1+m_2,7)/pow(m_1+m_2+m_3,3) * pow(m_3,7)/pow(m_1*m_2,3) * pow(L_1(m_1,m_2,a_1),4)/(pow(L_2(m_1,m_2,m_3,a_2),3) * pow(G_2(m_1,m_2,m_3,a_2,e_2),3));
}

double epsilon_m(double m_1,double m_2,double a_1,double a_2,double e_2){
  return (m_1-m_2)/(m_1+m_2) * a_1/a_2 * e_2/(1.0-e_2*e_2);
}

// Ojo, si es negativo? ver eq 24 de Naoz 2013
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

  ////////////////////////////////////////////////////////////
  // bulk properties
  
  
  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  
  double da_in_dt_orb;
  double da_in_dt_tid;
 
  da_in_dt_orb = 0.0;
  da_in_dt_tid = -2.0*a_in*( W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)
			     + W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) - 
    2.0*a_in*pow(e_in,2)/(1.0-e_in*e_in) * ( V_A(tv_A,R_A,m_A,m_B,k_A,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + 
					     V_B(tv_B,R_B,m_A,m_B,k_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) );
  
  return (params.q_orb*da_in_dt_orb) + (params.q_tid*da_in_dt_tid);

}






double da_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties
 
  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  
  double da_out_dt_orb;
  double da_out_dt_tid;
 
  da_out_dt_orb = 0.0;
  da_out_dt_tid = 0.0;
  
  return (params.q_orb*da_out_dt_orb) + (params.q_tid*da_out_dt_tid);
}



double de_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double de_in_dt_orb;
  double de_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  de_in_dt_orb = C2(m_A,m_B,m_C,a_in,a_out,e_out)*(1-e_in*e_in)/G_1(m_A,m_B,a_in,e_in) * 30.0*e_in*pow(sin(I_tot),2) * sin(2.0*w_in) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_out*(1-e_in*e_in)/G_1(m_A,m_B,a_in,e_in)*(35.0*cphi*pow(sin(I_tot),2)*e_in*e_in*sin(2.0*w_in) -
										 10.0*cos(I_tot)*pow(sin(I_tot),2)*cos(w_in)*sin(w_out)*(1.0-e_in*e_in) -
										 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)));
  
  de_in_dt_tid = -(V_A(tv_A,R_A,m_A,m_B,k_A,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + V_B(tv_B,R_B,m_A,m_B,k_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz))*e_in;
  
  return (params.q_orb*de_in_dt_orb) + (params.q_tid*de_in_dt_tid);
}



double de_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double de_out_dt_orb;
  double de_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  de_out_dt_orb = -C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*(1.0-e_out*e_out)/G_2(m_A,m_B,m_C,a_out,e_out) *
    (10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out) +
     A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)));
  de_out_dt_tid = 0.0;
  
  return (params.q_orb*de_out_dt_orb) + (params.q_tid*de_out_dt_tid);
}




double dI_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dI_in_dt_orb;
  double dI_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;
  double dGin_dt;
  double dGout_dt;
  double dHin_dt;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  

  dGin_dt  = -C2(m_A,m_B,m_C,a_in,a_out,e_out)*30.0*e_in*e_in*sin(2.0*w_in)*pow(sin(I_tot),2) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(-35.0*e_in*e_in*pow(sin(I_tot),2)*sin(2.0*w_in)*cphi +
						 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)) +
						 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*cos(w_in)*sin(w_out));
  
  dGout_dt = C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)) +
							  10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out));
  dHin_dt  = (sin(I_out)*dGin_dt - sin(I_in)*dGout_dt)/sin(I_tot);

  dI_in_dt_orb = -1.0 * ( dHin_dt - dGin_dt * cos(I_in) ) / (sin(I_in) * G_1(m_A,m_B,a_in,e_in) );
  dI_in_dt_tid = ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * cos(w_in) - 
                 ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * sin(w_in);

  
  return (params.q_orb*dI_in_dt_orb) + (params.q_tid*dI_in_dt_tid);
}



double dI_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dI_out_dt_orb;
  double dI_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;
  double dGin_dt;
  double dGout_dt;
  double dHin_dt;
  double dHout_dt;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dGin_dt  = -C2(m_A,m_B,m_C,a_in,a_out,e_out)*30.0*e_in*e_in*sin(2.0*w_in)*pow(sin(I_tot),2) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(-35.0*e_in*e_in*pow(sin(I_tot),2)*sin(2.0*w_in)*cphi +
						 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)) +
						 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*cos(w_in)*sin(w_out));
  dGout_dt = C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)) +
							  10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out));
  dHin_dt  = (sin(I_out)*dGin_dt - sin(I_in)*dGout_dt)/sin(I_tot);
  dHout_dt = -1.0*dHin_dt;
    
  dI_out_dt_orb = -1.0 * (dHout_dt - dGout_dt * cos(I_out))/(sin(I_out)*G_2(m_A,m_B,m_C,a_out,e_out));
  dI_out_dt_tid = 0.0;
  
  return (params.q_orb*dI_out_dt_orb) + (params.q_tid*dI_out_dt_tid);
}




double dW_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dW_in_dt_orb;
  double dW_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  dW_in_dt_orb = -3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in)) * (2.0+3.0*e_in*e_in-5.0*e_in*e_in*cos(2.0*w_in))*sin(2.0*I_tot) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(5.0*B*cos(I_tot)*cphi - A*sin(w_in)*sin(w_out) +
						 10.0*(1.0-3.0*pow(cos(I_tot),2))*(1.0-e_in*e_in)*sin(w_in)*sin(w_out))*sin(I_tot)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in));

  dW_in_dt_tid = ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) )
    * sin(w_in)/sin(I_in+I_out) +				
    ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) )
    * cos(w_in)/sin(I_in+I_out);

  return (params.q_orb * dW_in_dt_orb) + (params.q_tid * dW_in_dt_tid);

}




// Ojo aqui

double dW_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dW_out_dt_orb;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dW_out_dt_orb = -3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in)) * (2.0+3.0*e_in*e_in-5.0*e_in*e_in*cos(2.0*w_in))*sin(2.0*I_tot) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(5.0*B*cos(I_tot)*cphi - A*sin(w_in)*sin(w_out) +
						 10.0*(1.0-3.0*pow(cos(I_tot),2))*(1.0-e_in*e_in)*sin(w_in)*sin(w_out))*sin(I_tot)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in));

  return (params.q_orb * dW_out_dt_orb);
  
}


double dw_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  double cvel = 300.0e6 * params.uT/params.uL;
  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dw_in_dt_orb;
  double dw_in_dt_tid;
  double dw_in_dt_GR;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dw_in_dt_orb = 6.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)*(1.0/G_1(m_A,m_B,a_in,e_in)*(4.0*pow(cos(I_tot),2)
										   + (5.0*cos(2.0*w_in)-1.0)*(1.0-e_in*e_in-pow(cos(I_tot),2))) + 
						       cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out) * (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in)))) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_out*(e_in*(1.0/G_2(m_A,m_B,m_C,a_out,e_out) + cos(I_tot)/G_1(m_A,m_B,a_in,e_in)) *
					    (sin(w_in)*sin(w_out)*(10.0*(3.0*pow(cos(I_tot),2)-1.0)*(1.0-e_in*e_in)+A) - 5.0*B*cos(I_tot)*cphi) -
					    (1.0-e_in*e_in)/(e_in*G_1(m_A,m_B,a_in,e_in))*(sin(w_in)*sin(w_out)*10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-3.0*e_in*e_in) +
											   cphi*(3.0*A-10.0*pow(cos(I_tot),2)+2.0))); 

  
  dw_in_dt_tid = Z_A(R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Z_B(R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)  - 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * sin(w_in)*cos(I_in+I_out)/sin(I_in+I_out) - 
    ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * cos(w_in)*cos(I_in+I_out)/sin(I_in+I_out);
		

  dw_in_dt_GR = 3.0*pow(G*(m_A+m_B),1.5)/(pow(a_in,2.5)*cvel*cvel*(1.0-e_in*e_in));

  return (params.q_orb * dw_in_dt_orb) + (params.q_tid * dw_in_dt_tid) + (params.q_GR * dw_in_dt_GR);


}




double dw_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  
  double dw_out_dt_orb;
  double dw_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dw_out_dt_orb = 3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)*(2.0*cos(I_tot)/G_1(m_A,m_B,a_in,e_in) * (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in))) + 
							1.0/G_2(m_A,m_B,m_C,a_out,e_out) * (4.0+6.0*e_in*e_in+(5.0*pow(cos(I_tot),2)-3.0)*
											    (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in))))) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*(sin(w_in)*sin(w_out)*((4.0*e_out*e_out + 1.0)/(e_out*G_2(m_A,m_B,m_C,a_out,e_out))*
								 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in) -
								 e_out*(1.0/G_1(m_A,m_B,a_in,e_in) + cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out)) *
								 (A+10.0*(3.0*pow(cos(I_tot),2)-1.0)*(1.0-e_in*e_in))) +
					   cphi*(5.0*B*cos(I_tot)*e_out*(1.0/G_1(m_A,m_B,a_in,e_in) + cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out)) +
						 (4.0*e_out*e_out + 1.0)/(e_out*G_2(m_A,m_B,m_C,a_out,e_out))*A));
  
  dw_out_dt_tid = 0.0;
  
  return (params.q_orb * dw_out_dt_orb) + (params.q_tid * dw_out_dt_tid);

}





double dOm_Ax_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_A = params.gyr_rad_A;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid*( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-cos(W_in)*sin(w_in) - sin(W_in)*cos(w_in)*cos(I_in)) + 
      W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*sin(W_in)) - 
      Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(cos(W_in)*cos(w_in)-sin(W_in)*sin(w_in)*cos(I_in))));

}




double dOm_Ay_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_A = params.gyr_rad_A;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-sin(W_in)*sin(w_in) + cos(W_in)*cos(w_in)*cos(I_in)) + 
      W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-sin(I_in)*cos(W_in)) - 
      Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(W_in)*cos(w_in)+cos(W_in)*sin(w_in)*cos(I_in))) ) ;

}





double dOm_Az_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_A = params.gyr_rad_A;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);
  
  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
			  ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*cos(w_in)) + 
			    W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*cos(I_in) - 
			    Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*sin(w_in))) );

}





double dOm_Bx_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_B = params.gyr_rad_B;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid*( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
    ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-cos(W_in)*sin(w_in) - sin(W_in)*cos(w_in)*cos(I_in)) + 
      W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*sin(W_in)) - 
      Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(cos(W_in)*cos(w_in)-sin(W_in)*sin(w_in)*cos(I_in))));

}




double dOm_By_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_B = params.gyr_rad_B;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
    ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-sin(W_in)*sin(w_in) + cos(W_in)*cos(w_in)*cos(I_in)) + 
      W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-sin(I_in)*cos(W_in)) - 
      Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(W_in)*cos(w_in)+cos(W_in)*sin(w_in)*cos(I_in))) ) ;

}





double dOm_Bz_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
  ////////////////////////////////////////////////////////////
  // bulk properties

  
  double tv_A = params.tv_A;
  double tv_B = params.tv_B;
  double R_A  = params.R_A;
  double R_B  = params.R_B;
  double k_A  = params.k_A;
  double k_B  = params.k_B;
  double m_A  = params.m_A;
  double m_B  = params.m_B;
  double m_C  = params.m_C;
  double gyr_rad_B = params.gyr_rad_B;
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);
  
  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
			  ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*cos(w_in)) + 
			    W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*cos(I_in) - 
			    Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*sin(w_in))) );

}




////////////////////////////////////////////////////////////
//
// State vector function
//
////////////////////////////////////////////////////////////


int func (double t, const double y[], double f[], void *params){
  (void)(t); // avoid unused parameter warning

  Inpar parameters = *(Inpar *)params;
    
  double a_in  = y[0];
  double a_out = y[1];
  double e_in  = y[2];
  double e_out = y[3];
  double I_in  = y[4];
  double I_out = y[5];
  double W_in  = y[6];
  double W_out = y[7];
  double w_in  = y[8];
  double w_out = y[9];
  double Om_Ax = y[10];
  double Om_Ay = y[11];
  double Om_Az = y[12];
  double Om_Bx = y[13];
  double Om_By = y[14];
  double Om_Bz = y[15];

  f[0]  = da_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[1]  = da_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[2]  = de_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[3]  = de_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[4]  = dI_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[5]  = dI_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[6]  = dW_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[7]  = dW_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[8]  = dw_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[9]  = dw_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[10] = dOm_Ax_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[11] = dOm_Ay_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[12] = dOm_Az_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[13] = dOm_Bx_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[14] = dOm_By_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[15] = dOm_Bz_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  
  return GSL_SUCCESS;
  
}



/*
typedef struct my_estruct{
  double mu;
} exx;


double pos(double x, double vx, exx params){
  double mu = params.mu;
  return -mu*mu*x;
}

double vel(double x, double vx, exx params){
  return vx;
}

int func (double t, const double y[], double f[], void *params){
  (void)(t); // avoid unused parameter warning 
  //double mu = *(double *)params;
 
  exx stt   = *(exx *)params;
  stt.mu = 1.0;
  double x  = y[0];
  double vx = y[1];

 
  //f[0] = vx;
  //f[1] = -mu*mu*x;
 

  f[0] = vel(x,vx,stt);
  f[1] = pos(x,vx,stt);
  
  return GSL_SUCCESS;
}
*/

/*
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
*/
