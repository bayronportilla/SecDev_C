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
//#include "ModQuad.h"
//#include "matrix.h"
#include "ModOct.h"
#include "ConverCan.h"
#include <libconfig.h>

#define PI 3.14159


int main (void){

  FILE *fp;
  fp = fopen("data.dat","w");

  Inpar st;
  Inpar stc;
  st = params();
  stc = ConverToCan(st);

  //gsl_odeiv2_driver_alloc_y_new(const gsl_odeiv2_system * sys,
  //                              const gsl_odeiv2_step_type * T,
  //                              const double hstart,
  //                              const double epsabs,
  //                              const double epsrel)
  
  //const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkck;

  gsl_odeiv2_system sys = { func, NULL, 16, &stc}; // Define sistema de ecuaciones
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,1e-3, 1e-8, 1e-8);
    
  //const double h = 6.283e1;
  
  double y[16] = { stc.a_in, stc.a_out, stc.e_in, stc.e_out, 
		   stc.I_in, stc.I_out, stc.W_in, stc.W_out,
		   stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
		   stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz}; // y[number of entries of the array] = {}
  
   int i, s;
   double t = stc.t_ini;


   double progress;
   /*
   printf("m_A = %1.9e \n",stc.m_A);
   printf("m_B = %1.9e \n",stc.m_B);
   printf("m_C = %1.9e \n",stc.m_C);
   printf("R_A = %1.9e \n",stc.R_A);
   printf("R_B = %1.9e \n",stc.R_B);
   printf("P_rot_A = %1.9e \n",stc.P_rot_A);
   printf("P_rot_B = %1.9e \n",stc.P_rot_B);
   printf("alpha_A = %1.9e \n",stc.alpha_A);
   printf("alpha_B = %1.9e \n",stc.alpha_B);
   printf("beta_A = %1.9e \n",stc.beta_A);
   printf("beta_B = %1.9e \n",stc.beta_B);
   printf("tv_A = %1.9e \n",stc.tv_A);
   printf("tv_B = %1.9e \n",stc.tv_B);
   printf("k_A = %1.9e \n",stc.k_A);
   printf("k_B = %1.9e \n",stc.k_B);
   printf("a_in = %1.9e \n",stc.a_in);
   printf("a_out = %1.9e \n",stc.a_out);
   printf("e_in = %1.9e \n",stc.e_in);
   printf("e_out = %1.9e \n",stc.e_out);
   printf("I_in = %1.17e \n",stc.I_in);
   printf("I_out = %1.9e \n",stc.I_out);
   printf("W_in = %1.9e \n",stc.W_in);
   printf("W_out = %1.9e \n",stc.W_out);
   printf("w_in = %1.9e \n",stc.w_in);
   printf("w_out = %1.9e \n",stc.w_out);
   printf("Om_Ax = %1.9e \n",stc.Om_Ax);
   printf("Om_Ay = %1.9e \n",stc.Om_Ay);
   printf("Om_Az = %1.9e \n",stc.Om_Az);
   printf("Om_Bx = %1.9e \n",stc.Om_Bx);
   printf("Om_By = %1.9e \n",stc.Om_By);
   printf("Om_Bz = %1.9e \n",stc.Om_Bz);

   printf("%e\n",X_A(stc.tv_A,stc.R_A,stc.k_A,stc.m_A,
		     stc.m_B, stc.a_in, stc.e_in, stc.W_in,
		     stc.I_in, stc.w_in, stc.Om_Ax, stc.Om_Ay,
		     stc.Om_Az));
   
   */

   
   
   /*   
   printf("%e\n",Y_A(stc.tv_A,stc.R_A,stc.k_A,stc.m_A,
		     stc.m_B, stc.a_in, stc.e_in, stc.W_in,
		     stc.I_in, stc.w_in, stc.Om_Ax, stc.Om_Ay,
		     stc.Om_Az));
   */
   
   //double Om_A_in[3] = {stc.Om_Ax,stc.Om_Ay,stc.Om_Az};
   //printf("%e\n",RotMatt(stc.W_in,stc.I_in,stc.w_in,Om_A_in,"InOrb")[1]);

   /*
   double Om_A_in[3] = {stc.Om_Ax,stc.Om_Ay,stc.Om_Az};
   double Om_Ax_orb = RotMatt(stc.W_in,stc.I_in,stc.w_in,Om_A_in,"InOrb")[0];
   double Om_Ay_orb = RotMatt(stc.W_in,stc.I_in,stc.w_in,Om_A_in,"InOrb")[1];
   double Om_Az_orb = RotMatt(stc.W_in,stc.I_in,stc.w_in,Om_A_in,"InOrb")[2];

   printf("angulos: \n");
   printf("W=%e\n",stc.W_in);
   printf("I=%e\n",stc.I_in);
   printf("w=%e\n",stc.w_in);
   printf("\n");
   printf("valores iniciales: \n");
   printf("Om_Ax=%e\n",stc.Om_Ax);
   printf("Om_Ay=%e\n",stc.Om_Ay);
   printf("Om_Az=%e\n",stc.Om_Az);
   printf("\n");
   printf("valores finales: \n");
   printf("Om_Ax=%e\n",Om_Ax_orb);
   printf("Om_Ay=%e\n",Om_Ay_orb);
   printf("Om_Az=%e\n",Om_Az_orb);
   */

   /*
   printf("%e\n",tf_A(stc.tv_B,stc.R_B,stc.k_B,
		      stc.a_in,stc.m_A,stc.m_B));
   */
   
   //   W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az);

   /*
   printf("dI_in_dt = %e\n",dI_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				      stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				      stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				      stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				      0.0,stc));
   */
   printf("da_in_dt = %e\n",da_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				     stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				     stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				     stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				     0.0,stc));

   printf("da_out_dt = %e\n",da_out_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				      stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				      stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				      stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				      0.0,stc));

   printf("de_in_dt = %e\n",de_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				     stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				     stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				     stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				     0.0,stc));

   printf("de_out_dt = %e\n",de_out_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));

   printf("dI_in_dt = %e\n",dI_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				     stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				     stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				     stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				     0.0,stc));

   printf("dI_out_dt = %e\n",dI_out_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));

   printf("dW_in_dt = %e\n",dW_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				     stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				     stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				     stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				     0.0,stc));
   
   printf("dW_out_dt = %e\n",dW_out_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));

   printf("dw_in_dt = %e\n",dw_in_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				     stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				     stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				     stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				     0.0,stc));
   
   printf("dw_out_dt = %e\n",dw_out_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   printf("dOm_Ax_dt = %e\n",dOm_Ax_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   printf("dOm_Ay_dt = %e\n",dOm_Ay_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   printf("dOm_Az_dt = %e\n",dOm_Az_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));

   printf("dOm_Bx_dt = %e\n",dOm_Bx_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   printf("dOm_By_dt = %e\n",dOm_By_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   printf("dOm_Bz_dt = %e\n",dOm_Bz_dt(stc.a_in,stc.a_out,stc.e_in,stc.e_out,
				       stc.I_in, stc.I_out, stc.W_in, stc.W_out,
				       stc.w_in, stc.w_out, stc.Om_Ax, stc.Om_Ay,
				       stc.Om_Az, stc.Om_Bx, stc.Om_By, stc.Om_Bz,
				       0.0,stc));
   
   

   //printf("%f\n",st.gyr_rad_A);
   // exit(0);
   double ti = 1e3;
   while(t<stc.t_end){
     //s = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);
    
     s = gsl_odeiv2_driver_apply(d, &t, ti, y);
     if (s != GSL_SUCCESS){
       printf ("error: driver returned %d\n", s);
       break;
     }

     ti += 1e3;
     
    /*
    progress = (t/stc.t_end)*100.0;
    printf("Progress: %d per cent \n",(int)progress);
    
    if ( (int)progress%10 == 0){
      printf("Progress: %d per cent \n",(int)progress);
    }
    */
     printf("%.5e %.5e %.5e\n", t, y[0], y[2]);
     fprintf(fp,"%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n",
	     t,y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12],y[13],y[14],y[15]);
   }
   
   gsl_odeiv2_driver_free (d);
   fclose(fp);
   
   return 0;
   
}

