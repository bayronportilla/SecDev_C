#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h>

/* This example reads the configuration file 'example.cfg' and displays
 * some of its contents.
 */

//int main(int argc, char **argv){
int main(void){

  config_t cfg;
  config_setting_t *setting;
  const char *str;
  const char *met;
  int N;
  
  config_init(&cfg);

  // Read the file. If there is an error, report it and exit. 
  if(! config_read_file(&cfg, "example.cfg"))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    return(EXIT_FAILURE);
  }

  
  
  config_lookup_int(&cfg, "N", &N);
  printf("\nN = %d\n", N);

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
  double mass_A;
  double radius_A;
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

      config_setting_lookup_string(object, "name", &name_A);
      config_setting_lookup_float(object, "mass", &mass_A);
      config_setting_lookup_float(object, "radius", &radius_A);
      config_setting_lookup_float(object, "k", &k_A);
      config_setting_lookup_float(object, "tv", &tv_A);
      config_setting_lookup_float(object, "gyr_rad", &gyr_rad_A);
      config_setting_lookup_float(object, "a", &a_in);
      config_setting_lookup_float(object, "e", &e_in);
      config_setting_lookup_float(object, "W", &W_in);
      config_setting_lookup_float(object, "w", &w_in);
     
    }
    putchar('\n');
  }



  ////////////////////////////////////////////////////////////
  //
  // Getting information about inner B
  //
  ////////////////////////////////////////////////////////////
  
  const char *name_B;
  double mass_B;
  double radius_B;
  double k_B;
  double tv_B;
  double gyr_rad_B;

  
  setting = config_lookup(&cfg, "bodies.inner_B");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;

    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);

      config_setting_lookup_string(object, "name", &name_B);
      config_setting_lookup_float(object, "mass", &mass_B);
      config_setting_lookup_float(object, "radius", &radius_B);
      config_setting_lookup_float(object, "k", &k_B);
      config_setting_lookup_float(object, "tv", &tv_B);
      config_setting_lookup_float(object, "gyr_rad", &gyr_rad_B);

    }
}


  ////////////////////////////////////////////////////////////
  //
  // Getting information about outer body
  //
  ////////////////////////////////////////////////////////////
  
  const char *name_C;
  double mass_C;
  double radius_C;
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

      config_setting_lookup_string(object, "name", &name_C);
      config_setting_lookup_float(object, "mass", &mass_C);
      config_setting_lookup_float(object, "radius", &radius_C);
      config_setting_lookup_float(object, "k", &k_C);
      config_setting_lookup_float(object, "tv", &tv_C);
      config_setting_lookup_float(object, "gyr_rad", &gyr_rad_C);
      config_setting_lookup_float(object, "a", &a_out);
      config_setting_lookup_float(object, "e", &e_out);
      config_setting_lookup_float(object, "W", &W_out);
      config_setting_lookup_float(object, "w", &w_out);

    }
}

  config_destroy(&cfg);
  return(EXIT_SUCCESS);
}
