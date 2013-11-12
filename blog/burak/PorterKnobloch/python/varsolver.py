#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by velocityvar in twomode.py
"""

from scipy.integrate import odeint
from scipy.optimize import fsolve, newton
import twomode
import numpy as np

epsilon = np.finfo(np.float64).eps

#Callable function version:
def integrate(x0, p, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	"""
	xsol = odeint(twomode.velocityvar, x0, t, args=(p,), atol = abserror, rtol = relerror)
	
	return xsol

def funT(stoptime, x0, p):
	"""
	F(x) = f^t(x0) - x(0) = 0 for a periodic orbit
	"""
	abserr = 10*epsilon
	relerr = 100*epsilon
	numpoints = 1000
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

	xsol = odeint(twomode.vfullssp, x0, t, args=(p,), atol = abserr, rtol = relerr)
	Fx = np.linalg.norm(xsol[np.size(xsol,0)-1,:] - x0)
	#Fx = xsol[np.size(xsol,0)-1,:] - x0
	
	return Fx
	

if __name__ == "__main__":
	
	#Load parameters:
	p = np.loadtxt('data/parameters.dat')
	
	#Initial conditions:
	x10=-0.311102644920502
	y10=-0.311102644920502
	x20=-0.07020440068369531
	y20=-0.3862670595294275
	
	#Find the period
	
	x00 = [x10, y10, x20, y20]
	T = newton(funT, float(11.5485), args=(x00,p), tol=1000*epsilon)
	
	# ODE solver parameters
	abserr = 10*epsilon
	relerr = 100*epsilon
	stoptime = T
	numpoints = 5
	#stoptime = 11.549
	#numpoints = 500000

	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Pack up the initial conditions:
	x0=np.zeros(4+4*4)
	x0[0:4] = [x10,y10,x20,y20]
	x0[4:4+4*4] = np.identity(4).reshape(16) # Initial condition for the Jacobian is identity

	# Call the ODE solver
	xsol = odeint(twomode.velocityvar, x0, t, args=(p,), atol = abserr, rtol = relerr)
	
	
	#Print the solution
	np.savetxt('data/varsolution.dat', xsol)
	np.savetxt('data/vartime.dat', t)
