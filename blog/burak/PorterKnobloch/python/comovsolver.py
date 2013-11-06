#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vhatvphi on the symmetry
reduced manifold. Can be called as a script or imported as a module.
"""

import numpy as np 
from scipy.integrate import odeint
import twomode

T = twomode.generator()

def vperp(x,t,p):
	"""
    Perpendicular velocity function
    """
	tx = np.dot(T,x)
	vel = twomode.vfullssp(x,t,p) - (np.dot(twomode.vfullssp(x,t,p), tx)/np.dot(tx,tx))*tx
	return vel
	
abserr = 1.0e-8
relerr = 1.0e-6

#Callable function version:
def integrate(x, p, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	"""
	xsol = odeint(vperp, x, t, args=(p,), atol = abserror, rtol = relerror)
	
	return xsol


#Integrate only if the module is called as a script:

if __name__ == "__main__":
	
	#Load parameters:
	p = np.loadtxt('data/parameters.dat')

	#Initial conditions:

	x10=-1.32138709067
	y10=-1.32138709067
	x20= 0.26866738527 
	y20=-1.61705230412
	
	# ODE solver parameters
	abserr = 1.0e-8
	relerr = 1.0e-6
	stoptime = 500
	numpoints = 50000
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Pack up the initial conditions:
	x0 = [x10,y10,x20,y20]
	
	# Call the ODE solver
	xsol = odeint(vperp, x0, t, args=(p,), atol = abserr, rtol = relerr)
	
	#Print the solution
	for t1, x1 in zip(t,xsol):
		print t1, x1[0], x1[1], x1[2], x1[3]
