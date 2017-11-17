//#include <stdio.h>
//#include <stdlib.h>
#include <libconfig.h>


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
  double W_in;
  double w_in;

  // Inner B properties
  char   name_B[50];
  double m_B ;
  double R_B ;
  double k_B ;
  double tv_B ;
  double gyr_rad_B ;

  // Inner C properties
  char   name_C[50];
  double m_C;
  double R_C;
  double k_C;
  double tv_C;
  double gyr_rad_C;
  double a_out;
  double e_out;
  double W_out;
  double w_out;

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
  config_lookup_int(&cfg  , "q_orb"  , &q_orb);
  config_lookup_int(&cfg  , "q_tid"  , &q_tid);
  config_lookup_int(&cfg  , "q_GR"   , &q_GR);
  config_lookup_int(&cfg  , "rcu"    , &rcu);
  config_lookup_int(&cfg  , "rse"    , &rse);

  rest.n_steps = n_steps;
  rest.t_ini   = t_ini;
  rest.t_end   = t_end;
  rest.q_orb   = q_orb;
  rest.q_tid   = q_tid;
  rest.q_GR    = q_GR;

  
  //strcpy(rest.name,met);
    
  //printf("\nN = %d\n", N);
  
  /*
    if(config_lookup_string(&cfg, "method", &met))
    printf("method: %s\n\n", met);
    else
    fprintf(stderr, "No 'name' setting in configuration file.\n");
    
  */
  
  /*
    const char *name_files;
    const char *method;
    double t_ini,t_end,h_step;
    int q_orb,q_tid,q_GR;
    
    if(config_lookup_string(&cfg, "name_files", &name_files)
    && config_lookup_string(&cfg, "method", &method)
    && config_lookup_float(&cfg, "t_ini", &t_ini)
    && config_lookup_float(&cfg, "t_end", &t_end)
    && config_lookup_float(&cfg, "h_step", &h_step)
    && config_lookup_int(&cfg, "q_orb", &q_orb)
    && config_lookup_int(&cfg, "q_tid", &q_tid)
    && config_lookup_int(&cfg, "q_GR", &q_GR)){
    
    printf("%s \n",name_files);
    printf("%s \n",method);
    printf("%f \n",t_ini);
    printf("%d \n",q_orb);
    
    }
  */
  
  
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
      
    }
    //putchar('\n');
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
  
  config_destroy(&cfg);
  //return(EXIT_SUCCESS);

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
  rest.W_in = W_in;
  rest.w_in = w_in;

  rest.a_out = a_out;
  rest.e_out = e_out;
  rest.W_out = W_out;
  rest.w_out = w_out;
  
  
  
  return rest;

  
}
