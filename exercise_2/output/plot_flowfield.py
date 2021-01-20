#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 17})
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Mordern Roman']})

# define the parameters
Lx = 52.0
Ly = 52.0

# load the simulation data
x, y, vx, vy = np.loadtxt("./flowfield.dat", unpack=True, skiprows=0)

# set up the plot
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.quiver(x, y, vx, vy)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
if Lx == Ly: ax.set_aspect('equal')
ax.grid()
plt.tight_layout()
# plt.show()
plt.savefig("./flowfield.pdf")