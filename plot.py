#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

data      = np.loadtxt("data.dat")
time = data[:,0]
pos = data[:,1:2]
vel = data[:,2:3]

fig = plt.figure()
ax = plt.axes()
ax.plot(time,pos,'.')
ax.plot(time,vel,'+')
fig.savefig("pshp.png")
