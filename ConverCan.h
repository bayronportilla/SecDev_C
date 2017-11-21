

Inpar ConverToCan(Inpar st){

  ////////////////////////////////////////////////////////////
  // Defining constants

  double AU = 149.6e9;
  double MS = 1.989e30;
  double RS = 6.957e8;
  double MJ = 1.898e27;
  double RJ = 6.9911e7;
  double ME = 5.972e24;
  double RE = 6.371e6;
  double YEARS = 365.25*86400;
  double DAYS = 86400.0; 
  double InRad = 3.14159265359/180.0;
  
  
  //////////////////////////////////////////////////////////////////////
  // Finding canonical units
  
  double uM =  units(1.989e30,149.6e9,"uT")[0];
  double uL =  units(1.989e30,149.6e9,"uT")[1];
  double uT =  units(1.989e30,149.6e9,"uT")[2];

  //////////////////////////////////////////////////////////////////////
  // Converting initial structure values to SI

  double m_A   = st.m_A*MS;
  double m_B   = st.m_B*MS;
  double m_C   = st.m_C*MS;
  
  double R_A   = st.R_A*RS;
  double R_B   = st.R_B*RS;
  double R_C   = st.R_C*RS;
  
  double tv_A  = st.tv_A*YEARS;
  double tv_B  = st.tv_B*YEARS;
  double tv_C  = st.tv_C*YEARS;
  
  double a_in  = st.a_in*AU;
  double a_out = st.a_out*AU;

  double P_rot_A = st.P_rot_A*DAYS;
  double P_rot_B = st.P_rot_B*DAYS;
  
  double Om_Ax = st.Om_Ax/DAYS;
  double Om_Ay = st.Om_Ay/DAYS;
  double Om_Az = st.Om_Az/DAYS;
  
  double Om_Bx = st.Om_Bx/DAYS;
  double Om_By = st.Om_By/DAYS;
  double Om_Bz = st.Om_Bz/DAYS;
  
  double t_ini = st.t_ini*YEARS;
  double t_end = st.t_end*YEARS;
    
 
  ////////////////////////////////////////////////////////////
  // Converting initial structure values into canonical values

  Inpar canonical_st;

  canonical_st.t_ini     = t_ini/uT;
  canonical_st.t_end     = t_end/uT;
  canonical_st.q_orb     = st.q_orb;
  canonical_st.q_tid     = st.q_tid;
  canonical_st.q_GR      = st.q_GR;
  canonical_st.rcu       = st.rcu;
  canonical_st.rse       = st.rse;
 
  canonical_st.m_A       = m_A/uM;
  canonical_st.m_B       = m_B/uM;
  canonical_st.m_C       = m_C/uM;
    
  canonical_st.R_A       = R_A/uL;
  canonical_st.R_B       = R_B/uL;
  canonical_st.R_C       = R_C/uL;
  
  canonical_st.k_A       = st.k_A;
  canonical_st.k_B       = st.k_B;
  canonical_st.k_C       = st.k_C;
  
  canonical_st.tv_A      = tv_A/uT;
  canonical_st.tv_B      = tv_B/uT;
  canonical_st.tv_C      = tv_C/uT;
  
  canonical_st.gyr_rad_A = st.gyr_rad_A;
  canonical_st.gyr_rad_B = st.gyr_rad_B;
  canonical_st.gyr_rad_C = st.gyr_rad_C;
    
  canonical_st.a_in      = a_in/uL;
  canonical_st.e_in      = st.e_in;
  canonical_st.I_in      = st.I_in;
  canonical_st.W_in      = st.W_in;
  canonical_st.w_in      = st.w_in;
  canonical_st.P_rot_A   = st.P_rot_A/uT;
  canonical_st.P_rot_B   = st.P_rot_B/uT;
  canonical_st.alpha_A   = st.alpha_A;
  canonical_st.alpha_B   = st.alpha_B;
  canonical_st.beta_A    = st.beta_A;
  canonical_st.beta_B    = st.beta_B;
  

  canonical_st.a_out     = a_out/uL;
  canonical_st.e_out     = st.e_out;
  canonical_st.I_out     = st.I_out;
  canonical_st.W_out     = st.W_out;
  canonical_st.w_out     = st.w_out;

  canonical_st.I_tot     = st.I_tot;

  canonical_st.Om_Ax     = Om_Ax*uT;
  canonical_st.Om_Ay     = Om_Ay*uT;
  canonical_st.Om_Az     = Om_Az*uT;
  canonical_st.Om_Bx     = Om_Bx*uT;
  canonical_st.Om_By     = Om_By*uT;
  canonical_st.Om_Bz     = Om_Bz*uT;
    

  return canonical_st;
}

