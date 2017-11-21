//#include <stdio.h>
//#include <stdlib.h>
#include <libconfig.h>


double G_1(double m_1,double m_2,double a_1,double e_1);
double G_2(double m_1,double m_2,double m_3,double a_2,double e_2);
double *RotMat(double W, double I, double w, double *in_array, char dir[]);

typedef struct Inpar_params{

  // Inner A properties
  char   name_A[50];
  double m_A;
  double R_A;
  double k_A;
  double tv_A;
  double gyr_rad_A;
  double a_in;
  double e_in;
  double I_in;
  double W_in;
  double w_in;
  double P_rot_A;
  double alpha_A;
  double beta_A;

  // Inner B properties
  char   name_B[50];
  double m_B ;
  double R_B ;
  double k_B ;
  double tv_B ;
  double gyr_rad_B ;
  double P_rot_B;
  double alpha_B;
  double beta_B;

  // Inner C properties
  char   name_C[50];
  double m_C;
  double R_C;
  double k_C;
  double tv_C;
  double gyr_rad_C;
  double a_out;
  double e_out;
  double I_out;
  double W_out;
  double w_out;

  // Mutual inclination
  double I_tot;

  // Rotational velocities
  double Om_Ax;
  double Om_Ay;
  double Om_Az;
  double Om_Bx;
  double Om_By;
  double Om_Bz;
  
  
  // General properties
  double t_ini;
  double t_end;
  int n_steps;
  int q_orb;
  int q_tid;
  int q_GR;
  int rcu;
  int rse;
  
  
} Inpar; // Defining a new datatype *Inpar*
         //which is a structure of the tyoe *Inpar_params*


Inpar params(){
//double params(){
  
  Inpar rest; //defining the REturnSTucture

  config_t cfg;
  config_setting_t *setting;
  const char *str;
  const char *met;

  config_init(&cfg);
  
  // Read the file. If there is an error, report it and exit. 
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      //return(EXIT_FAILURE);
    }

  ////////////////////////////////////////////////////////////
  //
  // Getting and storing General Properties
  //
  ////////////////////////////////////////////////////////////

  int N;
  double t_ini;
  double t_end;
  double I_tot;
  int n_steps;
  int q_orb;
  int q_tid;
  int q_GR;
  int rcu;
  int rse;
  
  
  //config_lookup_string(&cfg, "method", &met);
  config_lookup_int(&cfg  , "n_steps", &n_steps);
  config_lookup_float(&cfg, "t_ini"  , &t_ini);
  config_lookup_float(&cfg, "t_end"  , &t_end);
  config_lookup_float(&cfg, "I_tot"  , &I_tot);
  config_lookup_int(&cfg  , "q_orb"  , &q_orb);
  config_lookup_int(&cfg  , "q_tid"  , &q_tid);
  config_lookup_int(&cfg  , "q_GR"   , &q_GR);
  config_lookup_int(&cfg  , "rcu"    , &rcu);
  config_lookup_int(&cfg  , "rse"    , &rse);


  
  ////////////////////////////////////////////////////////////
  //
  // Getting information about inner A
  //
  ////////////////////////////////////////////////////////////
  
  const char *name_A;
  double m_A;
  double R_A;
  double k_A;
  double tv_A;
  double gyr_rad_A;
  double a_in;
  double e_in;
  double W_in;
  double w_in;
  double P_rot_A;
  double alpha_A;
  double beta_A;
  
  
  setting = config_lookup(&cfg, "bodies.inner_A");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_string(object, "name",    &name_A);
      config_setting_lookup_float(object,  "mass",    &m_A);
      config_setting_lookup_float(object,  "radius",  &R_A);
      config_setting_lookup_float(object,  "k",       &k_A);
      config_setting_lookup_float(object,  "tv",      &tv_A);
      config_setting_lookup_float(object,  "gyr_rad", &gyr_rad_A);
      config_setting_lookup_float(object,  "a",       &a_in);
      config_setting_lookup_float(object,  "e",       &e_in);
      config_setting_lookup_float(object,  "W",       &W_in);
      config_setting_lookup_float(object,  "w",       &w_in);
      config_setting_lookup_float(object,  "P_rot_A", &P_rot_A);
      config_setting_lookup_float(object,  "alpha_A", &alpha_A);
      config_setting_lookup_float(object,  "beta_A",  &beta_A); 
    }
  }
  
  
  
  ////////////////////////////////////////////////////////////
  //
  // Getting information about inner B
  //
  ////////////////////////////////////////////////////////////
  
  const char *name_B;
  double m_B;
  double R_B;
  double k_B;
  double tv_B;
  double gyr_rad_B;
  double P_rot_B;
  double alpha_B;
  double beta_B;
  
  
  setting = config_lookup(&cfg, "bodies.inner_B");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_string(object, "name",    &name_B);
      config_setting_lookup_float(object,  "mass",    &m_B);
      config_setting_lookup_float(object,  "radius",  &R_B);
      config_setting_lookup_float(object,  "k",       &k_B);
      config_setting_lookup_float(object,  "tv",      &tv_B);
      config_setting_lookup_float(object,  "gyr_rad", &gyr_rad_B);
      config_setting_lookup_float(object,  "P_rot_B", &P_rot_B);
      config_setting_lookup_float(object,  "alpha_B", &alpha_B);
      config_setting_lookup_float(object,  "beta_B",  &beta_B);    
    }
  }
  
  
  ////////////////////////////////////////////////////////////
  //
  // Getting information about outer body
  //
  ////////////////////////////////////////////////////////////
  
  const char *name_C;
  double m_C;
  double R_C;
  double k_C;
  double tv_C;
  double gyr_rad_C;
  double a_out;
  double e_out;
  double W_out;
  double w_out;
  
  setting = config_lookup(&cfg, "bodies.outer");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_string(object, "name",    &name_C);
      config_setting_lookup_float(object,  "mass",    &m_C);
      config_setting_lookup_float(object,  "radius",  &R_C);
      config_setting_lookup_float(object,  "k",       &k_C);
      config_setting_lookup_float(object,  "tv",      &tv_C);
      config_setting_lookup_float(object,  "gyr_rad", &gyr_rad_C);
      config_setting_lookup_float(object,  "a",       &a_out);
      config_setting_lookup_float(object,  "e",       &e_out);
      config_setting_lookup_float(object,  "W",       &W_out);
      config_setting_lookup_float(object,  "w",       &w_out);
    }
  }


  ////////////////////////////////////////////////////////////
  //
  // Converting to radians
  //
  ////////////////////////////////////////////////////////////


  double InRad = M_PI/180.0;

  I_tot   = I_tot*InRad;
  W_in    = W_in*InRad;
  W_out   = W_out*InRad;
  w_in    = w_in*InRad;
  w_out   = w_out*InRad;
  alpha_A = alpha_A*InRad;
  beta_A  = beta_A*InRad;
  alpha_B = alpha_B*InRad;
  beta_B  = beta_B*InRad;

  
  ////////////////////////////////////////////////////////////
  //
  // Deriving inner and outer inclinations
  //
  ////////////////////////////////////////////////////////////


  double G_tot = pow( pow(G_1(m_A,m_B,a_in,e_in),2) + pow(G_2(m_A,m_B,m_C,a_out,e_out),2)
		      + 2*G_1(m_A,m_B,a_in,e_in)*G_2(m_A,m_B,m_C,a_out,e_out)*cos(I_tot) , 0.5 );
  
  double I_in     = asin(G_2(m_A,m_B,m_C,a_out,e_out)/G_tot * sin(I_tot));
  double I_out    = asin(G_1(m_A,m_B,a_in,e_in)/G_tot * sin(I_tot));


  ////////////////////////////////////////////////////////////
  //
  // Deriving initial rotational velocities for A and B
  //
  ////////////////////////////////////////////////////////////

  double Om_Ax_orb   = (2.0*M_PI/P_rot_A)*sin(beta_A)*cos(alpha_A);
  double Om_Ay_orb   = (2.0*M_PI/P_rot_A)*sin(beta_A)*sin(alpha_A);
  double Om_Az_orb   = (2.0*M_PI/P_rot_A)*cos(beta_A);
    
  double Om_Bx_orb   = (2.0*M_PI/P_rot_B)*sin(beta_B)*cos(alpha_B);
  double Om_By_orb   = (2.0*M_PI/P_rot_B)*sin(beta_B)*sin(alpha_B);
  double Om_Bz_orb   = (2.0*M_PI/P_rot_B)*cos(beta_B);

  double vels_A[3]   = {Om_Ax_orb,Om_Ay_orb,Om_Az_orb};
  double vels_B[3]   = {Om_Bx_orb,Om_By_orb,Om_Bz_orb};

  double Om_Ax_in = RotMat(W_in,I_in,w_in,vels_A,"OrbIn")[0];
  double Om_Ay_in = RotMat(W_in,I_in,w_in,vels_A,"OrbIn")[1];
  double Om_Az_in = RotMat(W_in,I_in,w_in,vels_A,"OrbIn")[2];
  double Om_Bx_in = RotMat(W_in,I_in,w_in,vels_B,"OrbIn")[0];
  double Om_By_in = RotMat(W_in,I_in,w_in,vels_B,"OrbIn")[1];
  double Om_Bz_in = RotMat(W_in,I_in,w_in,vels_B,"OrbIn")[2];
    
 
  // config_destroy(&cfg);

  rest.n_steps = n_steps;
  rest.t_ini   = t_ini;
  rest.t_end   = t_end;
  rest.q_orb   = q_orb;
  rest.q_tid   = q_tid;
  rest.q_GR    = q_GR;
  rest.rcu     = rcu;
  rest.rse     = rse;

  rest.m_A = m_A;
  rest.m_B = m_B;
  rest.m_C = m_C;

  rest.R_A = R_A;
  rest.R_B = R_B;
  rest.R_C = R_C;

  rest.k_A = k_A;
  rest.k_B = k_B;
  rest.k_C = k_C;

  rest.tv_A = tv_A;
  rest.tv_B = tv_B;
  rest.tv_C = tv_C;

  rest.gyr_rad_A = gyr_rad_A;
  rest.gyr_rad_B = gyr_rad_B;
  rest.gyr_rad_C = gyr_rad_C;

  rest.a_in = a_in;
  rest.e_in = e_in;
  rest.I_in = I_in;
  rest.W_in = W_in;
  rest.w_in = w_in;
  rest.P_rot_A = P_rot_A;
  rest.P_rot_B = P_rot_B;
  rest.alpha_A = alpha_A;
  rest.alpha_B = alpha_B;
  rest.beta_A  = beta_A;
  rest.beta_B  = beta_B;
  
  rest.a_out = a_out;
  rest.e_out = e_out;
  rest.I_out = I_out;
  rest.W_out = W_out;
  rest.w_out = w_out;

  rest.I_tot   = I_tot;

  rest.Om_Ax = Om_Ax_in;
  rest.Om_Ay = Om_Ay_in;
  rest.Om_Az = Om_Az_in;
  rest.Om_Bx = Om_Bx_in;
  rest.Om_By = Om_By_in;
  rest.Om_Bz = Om_Bz_in;
  
   
  return rest;
  config_destroy(&cfg);

}
