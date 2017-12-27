
#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np



sim_name  = "fig04_Naoz2013"
data_file = sim_name+".dat"
info_file = sim_name+".log"
fig_file  = sim_name+".png"
data  = np.loadtxt(data_file)
info  = np.genfromtxt(info_file,dtype=None)


uT = float(info[34][2])


time     = data[:,0]
a_in     = data[:,1:2]
a_out    = data[:,2:3]
e_in     = data[:,3:4]
e_out    = data[:,4:5]
I_in     = data[:,5:6]
I_out    = data[:,6:7]
W_in     = data[:,7:8]
W_out    = data[:,8:9]
w_in     = data[:,9:10]
w_out    = data[:,10:11]
Om_Ax_in = data[:,11:12]
Om_Ay_in = data[:,12:13]
Om_Az_in = data[:,13:14]
Om_Bx_in = data[:,14:15]
Om_By_in = data[:,15:16]
Om_Bz_in = data[:,16:17]



f,((ax_1), (ax_2)) = plt.subplots(2, 1,figsize=(10,8),sharex=True)
ax_1.plot(time*uT/(365.25*86400)/1e6,(I_in+I_out)*180.0/np.pi,'-',linewidth=2)
#ax_2.plot(time*uT/(365.25*86400)/1e6,1-e_in,'-',linewidth=2)
ax_2.plot(time*uT/(365.25*86400)/1e6,1-e_in,'.')#e,linewidth=2)
#ax_2.axhline(y=0.0023,xmin=0,xmax=4)
#ax_2.axvline(x=0.55,ymin=1e-5,ymax=1,color='k')
#ax_2.set_yscale("log")
#ax_1.set_xlim(0.0,100.0)
#ax_2.set_xlim(0.0,100.0)
#ax_1.set_ylim(0.0,180.0)
#ax_2.set_ylim(0.001,1.0)
f.savefig(fig_file)
