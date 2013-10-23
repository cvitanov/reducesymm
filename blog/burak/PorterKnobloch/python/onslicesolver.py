#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vinvpol in twomode.py
"""

import numpy as np 
from scipy.integrate import odeint
import twomode

#Load parameters:
p = np.loadtxt('data/parameters.dat')

T = twomode.generator()
xhatp = np.array([1,1,0,0],float)
tp = np.dot(T, xhatp)

def vphi(x,t,p):
	"""
    Velocity function for the group parameter
    """
	vel = np.dot(twomode.vfullssp(x,t,p),tp)/np.dot(np.dot(T,x),tp)
	return vel

def vhatvphi(xphi,t,p):
	"""
    Velocity function within the slice
    """
	vel=np.zeros(5)
	x = xphi[0:4]
	phi = xphi[4]
	vel[0:4] = twomode.vfullssp(x,t,p) - vphi(x,t,p)*np.dot(T,x)
	vel[4] = vphi(x,t,p)
	return vel


#Initial conditions:
x10 = -0.0101015007898    
x20 = -0.0101015007898
y10 = 0.00989852391435
y20 = 0.00990045343336
phi0 = 3.14169214968

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 1000
numpoints = 100000

# Create the time samples for the output of the ODE solver:
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the initial conditions:
xphi0 = [x10,x20,y10,y20, phi0]

# Call the ODE solver
xsol = odeint(vhatvphi, xphi0, t, args=(p,), atol = abserr, rtol = relerr)

#Print the solution
for t1, x1 in zip(t,xsol):
	print t1, x1[0], x1[1], x1[2], x1[3], x1[4]
