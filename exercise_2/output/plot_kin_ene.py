#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 17})
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Mordern Roman']})

# load the simulation data
h, KE = np.loadtxt("kin_ene.dat", unpack=True, skiprows=1)

# set up the plot
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.plot(h, KE, "g-", label=r"Simulation data")
ax.set_xlabel(r"$h$")
ax.set_ylabel(r"$E_K$")
ax.grid()
ax.legend()
plt.tight_layout()
# plt.show()
plt.savefig("kin_ene.pdf")
