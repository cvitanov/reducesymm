#
#multishoot.py
#
"""
Multiple shooting solver to find relative periodic orbits
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton, fsolve, root
import numdifftools as nd

import twomode
import sspsolver
import varsolver

#Read Porter - Knobloch system parameters:

pars = np.loadtxt('data/parameters.dat')

#Load the RPO candidates:

group = np.loadtxt('data/group.dat')
itineraries = np.loadtxt('data/itineraries.dat', dtype="str")
position = np.loadtxt('data/position.dat')
tof0 = np.loadtxt('data/tof0.dat')
phirpo0 = np.loadtxt('data/phi0.dat')
xrpo0 = np.loadtxt('data/xrpo0.dat')

#Let's handle the first one:
#Compute normal vector to the Poincare section:
sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
sectp3D = np.array([sectp[0], sectp[2], sectp[3]],float)
vaux = np.array([0,0,1], float) #Auxiliary vector to decide the Poincare section direction
unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float) #Unstable direction
nhat3D = np.cross(vaux, unstabledir)
nhat = np.array([nhat3D[0],0,nhat3D[1],nhat3D[2]], float)

#The slice:
T = twomode.generator()
xhatp = np.array([1,0,0,0],float)
tp = np.dot(T, xhatp)

#Initial guess:
x0 = xrpo0[0,:]
T0 = tof0[0]
phi0 = -phirpo0[0]

xt0 = np.append(x0, tof0[0])
xtphi0 = np.append(xt0, -phirpo0[0])

#Plotting modules:
from pylab import figure, grid, hold, legend, savefig
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

fig = plt.figure() #Create a figure instance
ax = fig.gca(projection = "3d") #Create an axis instance


stoptime = T0
numpoints = int(np.floor(T0/0.01))

t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

xsol = sspsolver.integrate(x0, pars, t, abserror=1.0e-14, relerror=1.0e-12)

ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], lw=0.5)
ax.set_xlabel('$x$', fontsize=18)
ax.set_ylabel('$y$', fontsize=18)
ax.set_zlabel('$z$', fontsize=18)

print "xsol[0,:]"
print xsol[0,:]

print "xsol[numpoints-1, :]"
print xsol[numpoints-1, :]

print "gxsol[numpoints-1,:]"
print np.dot(twomode.LieElement(phi0), xsol[numpoints-1,:])

plt.tight_layout()
plt.show()
