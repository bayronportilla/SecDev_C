#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

data      = np.loadtxt("data.dat")
time  = data[:,0]
e_in  = data[:,1:2]
I_tot = data[:,2:3]



"""
fig = plt.figure()
ax = plt.axes()
ax.plot(time,pos,'.')
ax.plot(time,vel,'+')
fig.savefig("pshp.png")
"""


f,((ax_1), (ax_2)) = plt.subplots(2, 1,figsize=(10,8),sharex=True)
ax_1.plot(time,e_in,'-',linewidth=2)
ax_2.plot(time,I_tot,'-',linewidth=2)
f.savefig("test.png")
