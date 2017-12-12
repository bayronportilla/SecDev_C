
#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

data  = np.loadtxt("data.dat")


time = data[:,0]
a_in = data[:,1:2]
a_out = data[:,2:3]
e_in = data[:,3:4]
e_out = data[:,4:5]
I_in = data[:,5:6]
I_out = data[:,6:7]
W_in = data[:,7:8]
W_out = data[:,8:9]
w_in = data[:,9:10]
w_out = data[:,10:11]
Om_Ax_in = data[:,11:12]
Om_Ay_in = data[:,12:13]
Om_Az_in = data[:,13:14]
Om_Bx_in = data[:,14:15]
Om_By_in = data[:,15:16]
Om_Bz_in = data[:,16:17]





"""
fig = plt.figure()
ax = plt.axes()
ax.plot(time,pos,'.')
ax.plot(time,vel,'+')
fig.savefig("pshp.png")
"""

f,((ax_1), (ax_2)) = plt.subplots(2, 1,figsize=(10,8),sharex=True)
ax_1.plot(time,(I_in+I_out)*180.0/np.pi,'-',linewidth=2)
ax_2.plot(time,1-e_in,'-',linewidth=2)
ax_2.set_yscale("log")
#ax_2.set_ylim(1e-5,1)
f.savefig("test.png")
