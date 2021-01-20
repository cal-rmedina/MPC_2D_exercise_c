#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 17})
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Mordern Roman']})

# define the parameters
h     = 0.1                    # MPC-collisoin time
n     = 5.0                    # fluid paticle density
kBT   = 1.0                    # thermal energy = temperature
alpha = np.pi/2.0              # collision angle
grav  = 0.0005                 # gravitational acceleration
Ly    = 52.0                   # channel hight

# kinetiv part of the (dynamic) viscosity
def eta_kin(h, n, alpha, T):
	f = n*T*h*(n/(1.0 - np.cos(2.0*alpha))/(n - 1.0 + np.exp(-n)) - 0.5)
	return f

# collisional part of the (dynamic) viscosity
def eta_coll(h, n, alpha):
	f = (1.0 - np.cos(alpha))/12.0/h*(n - 1.0 + np.exp(-n))
	return f

# full (dynamic) viscosity
def eta(h, n, alpha, T):
	f = eta_kin(h, n, alpha, T) + eta_coll(h, n, alpha)
	return f

# theoretical maximum velocity in channel
v_max = grav*Ly*Ly/(8.0*eta(h, n, alpha, kBT)/n)

# theoretical flowprofile
def v_x(y):
	f = 4.0*v_max*(Ly - y)*y/Ly/Ly
	return f

# setup the plotrange for theoretical flowprofile
y = np.arange(0, Ly, 0.1)

# load the simulation data
Y, V_X = np.loadtxt("./flowprofile.dat", unpack=True, skiprows=1)

# set up the plot
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.plot(Y, V_X, "bo-", label=r"Simulation data")
ax.plot(y, v_x(y), "r-", label=r"Theoretical result")
ax.set_xlabel(r"$y$")
ax.set_ylabel(r"$v_x$")
ax.grid()
ax.legend()
plt.tight_layout()
# plt.show()
plt.savefig("./flowprofile.pdf")