#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vinvpol in twomode.py
"""

import numpy as np 
from scipy.integrate import odeint
import twomode

#Parameter values:
mu1 = 1 
a1 = 0.47 
b1 = -1 
c1 = 1 
mu2 = -1 
a2 = 0 
b2 = 0 
c2 = -1 
e2 = 0
#mu1 = -2.8 
#a1 = -1 
#b1 = 0 
#c1 = -7.75 
#mu2 = 1 
#a2 = -2.66 
#b2 = 0 
#c2 = 1 
#e2 = 1

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
x10 = -0.384473125552
x20 = -0.384473125552
y10 = -1.22916029154
y20 = -1.50010256725
phi0 = 5.15722909618

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 100
numpoints = 10000

# Create the time samples for the output of the ODE solver:
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
xphi0 = [x10,x20,y10,y20, phi0]

# Call the ODE solver
xsol = odeint(vhatvphi, xphi0, t, args=(p,), atol = abserr, rtol = relerr)

#Print the solution
for t1, x1 in zip(t,xsol):
	print t1, x1[0], x1[1], x1[2], x1[3], x1[4]
