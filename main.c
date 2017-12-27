#include "allvars.h"


int main (void){

  Inpar st;
  st = params();

  //gsl_odeiv2_driver_alloc_y_new(const gsl_odeiv2_system * sys,
  //                              const gsl_odeiv2_step_type * T,
  //                              const double hstart,
  //                              const double epsabs,
  //                              const double epsrel)

  
  ////////////////////////////////////////////////////////////
  // This is the only line must be modified by the user
  //
  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  //
  ////////////////////////////////////////////////////////////

  
  gsl_odeiv2_system sys = { func, NULL, 16, &st}; // Define sistema de ecuaciones
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, 1e-3, 1e-8, 1e-8);



  double y[16] = { st.a_in, st.a_out, st.e_in, st.e_out, 
		   st.I_in, st.I_out, st.W_in, st.W_out,
		   st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
		   st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz}; // y[number of entries of the array] = {}
  
   int i, s;
   double t = st.t_ini;
   double progress;

   
   FetchInfo(st);

   /*
   printf("m_A = %1.9e \n",st.m_A);
   printf("m_B = %1.9e \n",st.m_B);
   printf("m_C = %1.9e \n",st.m_C);
   printf("R_A = %1.9e \n",st.R_A);
   printf("R_B = %1.9e \n",st.R_B);
   printf("P_rot_A = %1.9e \n",st.P_rot_A);
   printf("P_rot_B = %1.9e \n",st.P_rot_B);
   printf("alpha_A = %1.9e \n",st.alpha_A);
   printf("alpha_B = %1.9e \n",st.alpha_B);
   printf("beta_A = %1.9e \n",st.beta_A);
   printf("beta_B = %1.9e \n",st.beta_B);
   printf("tv_A = %1.9e \n",st.tv_A);
   printf("tv_B = %1.9e \n",st.tv_B);
   printf("k_A = %1.9e \n",st.k_A);
   printf("k_B = %1.9e \n",st.k_B);
   printf("a_in = %1.9e \n",st.a_in);
   printf("a_out = %1.9e \n",st.a_out);
   printf("e_in = %1.9e \n",st.e_in);
   printf("e_out = %1.9e \n",st.e_out);
   //printf("I_in = %1.17e \n",st.I_in);
   //printf("I_out = %1.9e \n",st.I_out);
   printf("W_in = %1.9e \n",st.W_in);
   printf("W_out = %1.9e \n",st.W_out);
   printf("w_in = %1.9e \n",st.w_in);
   printf("w_out = %1.9e \n",st.w_out);
   printf("Om_Ax_in = %1.9e \n",st.Om_Ax);
   printf("Om_Ay_in = %1.9e \n",st.Om_Ay);
   printf("Om_Az_in = %1.9e \n",st.Om_Az);
   printf("Om_Bx_in = %1.9e \n",st.Om_Bx);
   printf("Om_By_in = %1.9e \n",st.Om_By);
   printf("Om_Bz_in = %1.9e \n",st.Om_Bz);

   exit(0);
   */
   /*
   printf("%e\n",X_A(st.tv_A,st.R_A,st.k_A,st.m_A,
		     st.m_B, st.a_in, st.e_in, st.W_in,
		     st.I_in, st.w_in, st.Om_Ax, st.Om_Ay,
		     st.Om_Az));
   */



   //   exit(0);
   
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
   printf("dw_in_dt = %e\n",dw_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				      st.I_in, st.I_out, st.W_in, st.W_out,
				      st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				      st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				      0.0,st));
   */



   printf("da_in_dt = %1.17e\n",da_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				     st.I_in, st.I_out, st.W_in, st.W_out,
				     st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				     st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				     0.0,st));

     
   printf("da_out_dt = %1.17e\n",da_out_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				      st.I_in, st.I_out, st.W_in, st.W_out,
				      st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				      st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				      0.0,st));

   printf("de_in_dt = %1.17e\n",de_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				     st.I_in, st.I_out, st.W_in, st.W_out,
				     st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				     st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				     0.0,st));

   printf("de_out_dt = %1.17e\n",de_out_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));

   printf("dI_in_dt = %1.17e\n",dI_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				     st.I_in, st.I_out, st.W_in, st.W_out,
				     st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				     st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				     0.0,st));

   printf("dI_out_dt = %1.17e\n",dI_out_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));

   printf("dW_in_dt = %1.17e\n",dW_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				     st.I_in, st.I_out, st.W_in, st.W_out,
				     st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				     st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				     0.0,st));
   
   printf("dW_out_dt = %1.17e\n",dW_out_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));

   printf("dw_in_dt = %1.17e\n",dw_in_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				     st.I_in, st.I_out, st.W_in, st.W_out,
				     st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				     st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				     0.0,st));
   
   printf("dw_out_dt = %1.17e\n",dw_out_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));
   
   printf("dOm_Ax_dt = %1.17e\n",dOm_Ax_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));
   
   printf("dOm_Ay_dt = %1.17e\n",dOm_Ay_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));
   
   printf("dOm_Az_dt = %1.17e\n",dOm_Az_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));

   printf("dOm_Bx_dt = %1.17e\n",dOm_Bx_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));
   
   printf("dOm_By_dt = %1.17e\n",dOm_By_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));
   
   printf("dOm_Bz_dt = %1.17e\n",dOm_Bz_dt(st.a_in,st.a_out,st.e_in,st.e_out,
				       st.I_in, st.I_out, st.W_in, st.W_out,
				       st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
				       st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz,
				       0.0,st));

   
   
   
   FILE *fp;
   char src[100];
   char dest[100];
   char name_files[100];
   strcpy(src,st.sim_name);
   strcpy(dest,".dat");
   strcpy(name_files,strcat(src,dest));
   fp = fopen(name_files,"w");
   double ti = t;



   while(t<st.t_end){

     s = gsl_odeiv2_driver_apply(d, &t, ti, y);
     if (s != GSL_SUCCESS){
       printf ("error: driver returned %d\n", s);
       break;
     }
     ti += st.h_output;     
     
     progress = (t/st.t_end)*100.0;
     printf("Progress: %d per cent \n",(int)progress);
     
     if ( (int)progress%10 == 0){
       printf("Progress: %d per cent \n",(int)progress);
     }
     

     fprintf(fp,"%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n",
	     t,y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12],y[13],y[14],y[15]);
   }
   
   gsl_odeiv2_driver_free (d);
   fclose(fp);

   
   return 0;
   
}

